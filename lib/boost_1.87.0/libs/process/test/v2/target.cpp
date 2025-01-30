#include <iostream>
#include <string>
#include <thread>
#include <boost/asio/steady_timer.hpp>
#include <boost/process/v2/environment.hpp>

extern char **environ;

#if defined(BOOST_PROCESS_V2_WINDOWS)
#include <windows.h>
#else
#include <unistd.h>
#endif


int main(int argc, char * argv[])
{
    std::string mode = argv[1];
    if (mode == "exit-code")
        return std::stoi(argv[2]);
    else if (mode == "sleep")
    {
        const auto delay = std::chrono::milliseconds(std::stoi(argv[2]));
        std::this_thread::sleep_for(delay);
        return 0;
    }
    else if (mode == "print-args")
        for (auto i = 0; i < argc; i++)
        {
            std::cout << argv[i] << std::endl;
            std::cerr << argv[i] << std::endl;
            if (!std::cout || !std::cerr)
                return 1;
        }
    else if (mode == "echo")
        std::cout << std::cin.rdbuf();
    else if (mode == "print-cwd")
    {
#if defined(BOOST_PROCESS_V2_WINDOWS)
        wchar_t buf[65535];
        const auto sz = ::GetCurrentDirectoryW(sizeof(buf), buf);
        std::wcout << boost::process::v2::wstring_view(buf, sz) << std::flush;
#else
        char buf[65535];
        printf(::getcwd(buf, sizeof(buf)));
#endif
        return 0;
    }
    else if (mode == "check-eof")
    {
        std::string st;
        std::cin >> st;
        return std::cin.eof() ? 0 : 1;
    }
    else if (mode == "print-env")
    {
        auto p = ::getenv(argv[2]);
        if (p && *p)
            printf("%s", p);
        else
        {
            printf("Can't find %s in environment\n", argv[2]);
            for (auto e = environ; e != nullptr; e++)
                printf("    %s\n", *e);
            return 3;
        }
    }
#if defined(BOOST_PROCESS_V2_WINDOWS)
    else if (mode == "showwindow")
    {
        STARTUPINFO si;
        GetStartupInfo(&si);   
        return static_cast<int>(si.wShowWindow);
    }
#endif
    else if (mode == "sigterm")
    {
      boost::asio::io_context ctx;
      boost::asio::steady_timer tim{ctx, std::chrono::seconds(10)};
      static boost::asio::steady_timer * tim_p = &tim;

#if defined(BOOST_PROCESS_V2_WINDOWS)
      SetConsoleCtrlHandler(
          [](DWORD kind)
          {
            if (kind == CTRL_CLOSE_EVENT)
            {
              // windows doesn't like us doing antyhing else
              ::exit(0);
              if (tim_p != nullptr)
                tim_p->cancel();
              return TRUE;
            }
            else
              return FALSE;
          }, TRUE);
#else
      signal(SIGTERM, [](int) { if (tim_p != nullptr) tim_p->cancel();});
#endif

      boost::system::error_code ec;
      tim.async_wait([&](boost::system::error_code ec_) { ec = ec_; });
      ctx.run();
      tim_p = nullptr;
      return ec ? EXIT_SUCCESS : 32;
    }
    else if (mode == "sigint")
    {
      boost::asio::io_context ctx;
      boost::asio::steady_timer tim{ctx, std::chrono::seconds(10)};
      static boost::asio::steady_timer * tim_p = &tim;

#if defined(BOOST_PROCESS_V2_WINDOWS)
      SetConsoleCtrlHandler(NULL, FALSE);
      auto res = SetConsoleCtrlHandler(
              [](DWORD kind)
              {
                if (kind == CTRL_C_EVENT)
                {
                  if (tim_p != nullptr)
                    tim_p->cancel();
                  return TRUE;
                }
                else
                  return FALSE;
              }, TRUE);
      BOOST_ASSERT(res != FALSE);
#else
      signal(SIGINT, [](int) { if (tim_p != nullptr) tim_p->cancel();});
#endif

      boost::system::error_code ec;
      tim.async_wait([&](boost::system::error_code ec_) { ec = ec_; });
      ctx.run();
      tim_p = nullptr;
      return ec ? EXIT_SUCCESS : 33;
    }
    else
        return 34;

    return 0;
}