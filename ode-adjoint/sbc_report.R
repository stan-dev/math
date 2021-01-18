#' ---
#' title: Simulation based calibration for CVODES solvers
#' author: ""
#' date: "`r date()`"
#' output: html_vignette
#' params:
#'   include_plots: FALSE
#' ---
#'
#+ include=FALSE
source("sbc_tools.R")
library(knitr)
library(tools)
library(assertthat)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
theme_set(theme_bw())

knitr::opts_chunk$set(
  fig.width = 1.62*4,
  fig.height = 4,
  cache=FALSE,
  echo=FALSE
    )

calibration <- readRDS("calibration.rds")

plot_binned <- function(cal_df) {

    pl <- NULL

    if(!include_plots)
        return(pl)
    
    if(!all(cal_df$count == 0)) {
        
        S <- calibration$S
        B <- calibration$B
        
        c95 <- qbinom(c(0.025, 0.5, 0.975), S, 1 / B)
        
        dd <- cal_df %>%
            arrange(param, bin) %>%
                group_by(param) %>%
                    mutate(ecdf = cumsum(count) / S, ecdf_ref = (bin + 1) / B) %>%
                        filter(!all(ecdf == 0))
        
        nparam <- length(unique(dd$param))
        if(unique(dd$partype) %in% c("mu_eta", "tau_eta")){
            nc <- nparam
        } else{
            nc <- 2
        }
        
        nr <- max(1, ceiling(nparam / nc))

        pl <- list()
        pl[["hist"]] <- ggplot(dd, aes(bin, count)) +
            facet_wrap(~ param, nrow = nr, ncol = nc) +
                geom_col() +
                    geom_hline(yintercept=c95[c(1,3)], linetype=I(2)) +
                        geom_hline(yintercept=c95[c(2)], linetype=I(3))
        pl[["ecdf_diff"]] <- ggplot(dd, aes(bin, ecdf-ecdf_ref)) +
            facet_wrap(~ param, nrow = nr, ncol = nc) +
                geom_step() +
                    geom_hline(yintercept=0, linetype=I(3))
        pl
    }

    return(pl)
}


dim(calibration)
rank_params <- names(calibration)[grepl(names(calibration), pattern = "rank")]

# calibration_data_binned <- calibration_data[, scale64(.SD), by=c("problem", "model", params)]

sampled_ranks <- 2^8
B <- sampled_ranks/ 2^3
B

calibration_binned <- calibration %>%
    mutate_at(.vars = rank_params, .funs = function(x) ceiling((x + 1) / (sampled_ranks/ B) - 1)) %>%
        mutate(adjoint_integrator=factor(adjoint_integrator))


long <- calibration_binned[c("system_size", "adjoint_integrator", rank_params)] %>%
    tidyr::gather(key = "param", value = "rank", -system_size, -adjoint_integrator) %>%
        group_by(system_size, adjoint_integrator, param)  %>%
            mutate(all_NA=all(is.na(rank))) %>%
                filter(!all_NA)

metrics_long <- calibration_binned[!(names(calibration_binned) %in% rank_params)] %>%
  tidyr::gather(key = "metric", value = "value", -system_size, -adjoint_integrator, -job.id, -problem, -algorithm) 

long_split <- split(long, long$system_size)

pl_rank <- lapply(long_split, function(d) ggplot(d, aes(rank)) + facet_grid(param~adjoint_integrator*system_size, labeller=label_both) + geom_bar() + theme(strip.text.y = element_text(angle = 0)))

pl_rank[[1]] + theme(strip.text.y = element_text(angle = 0))
pl_rank[[2]]

ggsave("rank_size_1.png", pl_rank[[1]], width=6, height=6)
ggsave("rank_size_4.png", pl_rank[[2]], width=6, height=14)

metrics_split <- split(metrics_long, metrics_long$system_size)

pl_metrics <- lapply(metrics_split, function(d) ggplot(d, aes(adjoint_integrator, value)) +
    facet_wrap(~ metric*system_size, labeller=label_both, scales="free_y") + geom_boxplot()+
        theme(strip.text.y = element_text(angle = 0)))

ggsave("metrics_1.png", pl_metrics[[1]], width=10, height=10)
ggsave("metrics_4.png", pl_metrics[[2]], width=10, height=10)
