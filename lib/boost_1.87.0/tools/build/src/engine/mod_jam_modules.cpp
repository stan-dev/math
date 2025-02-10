/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_jam_modules.h"

#include "bindjam.h"
#include "compile.h"
#include "frames.h"
#include "function.h"
#include "hash.h"
#include "pathsys.h"
#include "rules.h"
#include "startup.h"
#include "variable.h"

#include "mod_jam_errors.h"
#include "mod_path.h"

#include <cassert>

namespace b2 { namespace jam { namespace modules {

list_ref binding(std::string module_name)
{
    return *jam::variable("modules", module_name + ".__binding__");
}

void record_binding(std::string module_name, value_ref value)
{
    std::string loading = jam::variable("modules", ".loading")[-1]->str();
    jam::variable("modules", loading + ".__binding__") = value;
}

void poke(std::string module_name, list_cref variables, list_cref value)
{
    for (auto variable : variables)
    {
        jam::variable(module_name, variable->str()) = value;
    }
}

list_ref peek(std::string module_name, list_cref variables)
{
    list_ref result;
    for (auto variable : variables)
    {
        result.append(jam::variable(module_name, variable->str()));
    }
    return result;
}

static void clone_rule(rule_ptr source_rule, module_ptr target_module)
{
    // IMPORT $(source-module) : $(source-rule) : $(target-module) :
    // $(source-rule) ;
    rule_ptr target_rule
        = import_rule(source_rule, target_module, source_rule->name);
    rule_localize(target_rule, target_module);
    // EXPORT $(target-module) : $(target-rule) ;
    target_rule->exported = 1;
    // IMPORT $(target-module) : $(target-rule) : :
    // $(target-module).$(target-rule) ;
    std::string target_rule_name = value_ref(target_rule->name)->str();
    std::string target_module_name = value_ref(target_module->name)->str();
    std::string global_rule_name = target_module_name + "." + target_rule_name;
    import_rule(target_rule, root_module(), value_ref(global_rule_name));
}

void clone_rules(std::tuple<std::string, std::string> source_target_modules)
{
    module_ptr source_mod
        = bindmodule(value_ref(std::get<0>(source_target_modules)));
    module_ptr target_mod
        = bindmodule(value_ref(std::get<1>(source_target_modules)));
    hash_enumerate(source_mod->rules, clone_rule, target_mod);
}

list_ref call_in(value_ref module_name,
    std::tuple<value_ref, list_ref> rule_name_a1,
    const lists & rest,
    bind::context_ref_ context_ref)
{
    value_ref rule_name = std::get<0>(rule_name_a1);
    frame * outer = context_ref.get<jam_context>().frame;
    module_scope scope(outer, module_name->str());
    lists args;
    args | *std::get<1>(rule_name_a1) | rest;
    LIST * result = call_rule(rule_name, outer, args);
    return list_ref(result, true);
}

list_ref call_locally(std::tuple<std::string, list_ref> rule_name_a1,
    const lists & rest,
    bind::context_ref_ context_ref)
{
    std::string rule_name = std::get<0>(rule_name_a1);
    auto dot_index = rule_name.find('.');
    std::string module_name;
    if (dot_index != std::string::npos)
    {
        module_name = rule_name.substr(0, dot_index);
        rule_name = rule_name.substr(dot_index + 1);
    }
    frame * outer = context_ref.get<jam_context>().frame;
    module_scope scope(outer, module_name.c_str());
    lists args;
    args | *std::get<1>(rule_name_a1) | rest;
    LIST * result = call_rule(value_ref(rule_name), outer, args);
    return list_ref(result, true);
}

void run_tests(value_ref m, bind::context_ref_ context_ref)
{
    variable tested_modules("modules", ".tested");
    list_cref argv = variable("ARGV");
    bool quiet = argv.contains("--quiet");
    bool debug = argv.contains("--debug");
    bool debug_module = argv.contains("--debug-module=" + m);
    bool debug_tests = argv.contains("--debug-tests");

    // Avoid recursive test invocations.
    if ((*tested_modules).contains(m)) return;
    // Only test if asked for in CLI.
    if (!debug && !debug_module) return;
    // Track what we've tested to avoid recursion.
    tested_modules += m;
    // Warn if we have no __test__ rule to run.
    list_ref rules(module_rules(bindmodule(m)), true);
    if (!rules.contains("__test__"))
    {
        if (!quiet && debug_tests)
            out_printf(
                "warning: no __test__ rule defined in module %s\n", m->str());
        // out_printf("\t%s rules: ", m->str());
        // list_print(*rules);
        // out_printf("\n");
        return;
    }
    // Debug info.
    if (!quiet) out_printf("testing module %s...\n", m->str());
    // Create a shell __test-M__ modules to run the tests in.
    module_ptr source_module = bindmodule(m);
    std::string test_module_name = std::string("__test-") + m->str() + "__";
    module_ptr test_module = bindmodule(value_ref(test_module_name));
    for (auto rule : rules)
        import_rule(find_rule(rule, source_module), test_module, rule);
    rule_localize(import_rule(find_rule(value_ref("__test__"), source_module),
                      test_module, value_ref("__test__")),
        test_module);
    // We can now run the tests in the __test-M__.__test__ rule.
    module_scope test_scope(
        context_ref.get<jam_context>().frame, test_module_name.c_str());
    run_rule(test_scope.frame(), "__test__");
}

void load(value_ref module_name,
    value_ref filename,
    list_cref search,
    bind::context_ref_ context_ref)
{
    jam_context & context = context_ref.get<jam_context>();
    variable loaded_v("modules", ".loaded");
    variable loading_v("modules", ".loading");
    variable untested_v("modules", ".untested");
    variable tested_v("modules", ".tested");
    std::string module_basename = PATHNAME(module_name->str()).base();
    if (!(*loaded_v).contains(module_name))
    {
        // Not already loaded.
        filename = filename.has_value() ? filename : module_name + ".jam";
        // Mark the module loaded so we do not try to load it recursively.
        loaded_v += module_basename;
        // Suppress tests if any module loads are already in progress.
        bool suppress_test = loading_v;
        // Push this module on the loading stack.
        loading_v += module_name;
        // Remember that it is untested.
        untested_v += module_name;
        // Insert the new module's __name__ and __file__ globals.
        variable(module_name, "__name__") = module_name;
        variable(module_name, "__file__") = filename;
        // Use BOOST_BUILD_PATH env as default.
        list_ref module_target_search
            = search.empty() ? *variable("BOOST_BUILD_PATH") : search;
        // Include the file, if we find it, and with some setup.
        {
            {
                // Add some grist so that the module will have a unique target
                // name.
                value_ref module_target = "<module@>" + filename;
                target_ref module_target_r(module_target);
                // Set the search path.
                module_target_r.on_set("SEARCH", module_target_search);
                // Set the binding rule.
                module_target_r.on_set("BINDRULE", "modules.record-binding");
                // Do the load, in the module.
                {
                    module_scope_in_function in_module(
                        context.frame, module_name->str());
                    parse_include(module_target, in_module.frame());
                    // Allow the module to see its own names with full
                    // qualification.
                    list_ref to_import(module_rules(in_module), true);
                    for (auto rule : to_import)
                    {
                        import_rule(find_rule(rule, in_module), in_module,
                            module_name + "." + rule);
                    }
                }
            }
            // Did the include succeed?
            if (module_name != "modules" && binding(module_name).empty())
            {
                run_rule(context.frame, "import", list_ref() + "errors");
                run_rule(context.frame, "errors.error",
                    list_ref() + "Could not find module" + module_name + "in"
                        + module_target_search);
            }
            // Pop the loading stack. Must happen before testing or we will run
            // into a circular loading dependency.
            loading_v = loading_v.get().slice(0, -2);
            // Run any pending tests if this is an outer load.
            if (!suppress_test)
            {
                for (auto to_test : list_ref(*untested_v))
                {
                    run_tests(to_test, context_ref);
                }
                untested_v = list_ref();
            }
        }
    }
    else if (!(*tested_v).contains(module_name))
    {
        // Native modules are pre-loaded but can be untested. In that case we
        // make sure we run tests for those.
        untested_v += module_name;
        if (!loading_v)
        {
            for (auto to_test : list_ref(*untested_v))
            {
                run_tests(to_test, context_ref);
            }
            untested_v = list_ref();
        }
    }
    else if ((*loading_v).contains(module_name))
    {
        run_rule(context.frame, "import", list_ref() + "errors");
        run_rule(context.frame, "errors.error",
            list_ref() + "loading \"" + module_name->str() + "\n",
            list_ref() + "circular module loading dependency:",
            list_ref() + *loading_v + "==>" + module_name);
    }
}

namespace detail {
value_ref caller_module(FRAME * frame, int levels = 0)
{
    for (; levels >= 0; --levels) frame = frame->prev ? frame->prev : frame;
    return frame->module == root_module() ? nullptr : frame->module->name;
}
} // namespace detail

void import(list_cref module_names,
    list_cref rules_opt,
    list_cref rename_opt,
    bind::context_ref_ context_ref)
{
    jam_context & context = context_ref.get<jam_context>();
    value_ref cwd = variable("modules", ".cwd")[0];
    variable tested_v("modules", ".tested");
    variable untested_v("modules", ".untested");
    bool wildcard_rules_opt = rules_opt == (list_ref() + "*");
    if ((wildcard_rules_opt || rules_opt.empty()) && !rename_opt.empty())
    {
        run_rule(context.frame, "import", list_ref() + "errors");
        run_rule(context.frame, "errors.error",
            list_ref()
                + "Rule aliasing is only available for explicit imports.");
    }
    if (module_names.length() >= 2
        && (!rules_opt.empty() || !rename_opt.empty()))
    {
        run_rule(context.frame, "import", list_ref() + "errors");
        run_rule(context.frame, "errors.error",
            list_ref()
                + "When loading multiple modules, no specific rules or renaming is allowed.");
    }

    value_ref caller_module = detail::caller_module(context.frame);

    // Import each specified module
    for (value_ref m : module_names)
    {
        PATHNAME module_pathname(m->str());
        std::string module_basename = module_pathname.base();

        // If the importing module is not already in the BOOST_BUILD_PATH,
        // prepend it to the path.  We do not want to invert the search
        // order of modules that are already there.
        std::string caller_location;
        if (caller_module.has_value())
        {
            list_ref caller_binding(binding(caller_module->str()));
            if (!caller_binding.empty())
            {
                caller_location = paths::normalize(
                    PATHNAME(caller_binding[0]->str()).dir());
                caller_location = paths::normalize(
                    paths::rooted(caller_location, cwd->str()));
            }
        }

        list_ref search;
        for (auto p : *variable("BOOST_BUILD_PATH"))
            search.push_back(paths::rooted(value_ref(p), cwd));

        if (!caller_location.empty() && !search.contains(caller_location))
            search = std::move(list_ref(caller_location).append(search));

        if (!module_pathname.dir().empty())
        {
            list_ref xsearch;
            xsearch.push_back(caller_location + "/" + module_pathname.dir());
            for (value_ref s : search)
                xsearch.push_back(std::string(s) + "/" + module_pathname.dir());
            xsearch.append(search);
            search = std::move(xsearch);
        }

        list_cref loading_r = *variable("modules", ".loading");
        list_cref loaded_r = *variable("modules", ".loaded");
        if (!loading_r.contains(module_basename))
        {
            if (!loaded_r.contains(module_basename))
            {
                load(module_basename, value_ref(), list_cref(*search),
                    context_ref);
            }
            else if (!(*tested_v).contains(m))
            {
                // Native modules are pre-loaded but can be untested. In that
                // case we make sure we run tests for those. The load takes care
                // of doing the testing.
                untested_v += m;
                load(module_basename, value_ref(), list_cref(*search),
                    context_ref);
            }
        }
        else if (!loaded_r.contains(module_basename) && !loading_r.empty()
            && loading_r[loading_r.length() - 1]->str() == module_basename)
        {
            // Special case for back-compat.
            b2::jam::errors::warning(lists()
                    | *(list_ref() + "loading" + module_basename)
                    | list_ref("circular module loading dependency:")
                    | *(list_ref() + loading_r + "==>" + module_basename),
                context_ref);
            load(module_basename, value_ref(), list_cref(*search), context_ref);
        }

        import_module(*list_ref(module_basename),
            caller_module.has_value() ? bindmodule(caller_module)
                                      : root_module());

        if (!rules_opt.empty())
        {
            list_ref source_names;
            if (wildcard_rules_opt)
                source_names.reset(module_rules(bindmodule(m)));
            else
                source_names.append(rules_opt);
            list_cref target_names(
                rename_opt.empty() ? *source_names : *rename_opt);
            run_rule(context.frame, "IMPORT", list_ref(module_basename),
                source_names, list_ref() + caller_module, target_names);
        }
    }
}

/*
Jam init code is:

# Copyright 2003 Dave Abrahams
# Copyright 2003, 2005 Vladimir Prus
# Copyright 2022 René Ferdinand Rivera Morell
*/
const char * modules_module::init_code = R"jam(

# Cache workdir on load.
.cwd = [ PWD ] ;

# Essentially an include guard; ensures that no module is loaded multiple times.
.loaded ?= ;

# A list of modules currently being loaded for error reporting of circular
# dependencies.
.loading ?= ;

# A list of modules needing to be tested using their __test__ rule.
.untested ?= ;

# A list of modules which have been tested using their __test__ rule.
.tested ?= ;

# These rules need to be available in all modules to implement module loading
# itself and other fundamental operations.
local globalize = peek poke record-binding ;
IMPORT modules : $(globalize) : : modules.$(globalize) ;
IMPORT modules : import : : import ;

rule __test__ ( )
{
    import assert ;

    module modules.__test__
    {
        foo = bar ;
    }

    assert.result bar : peek modules.__test__ : foo ;

    poke modules.__test__ : foo : bar baz ;
    assert.result bar baz : peek modules.__test__ : foo ;
}

)jam";

}}} // namespace b2::jam::modules
