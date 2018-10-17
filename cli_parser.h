#ifndef CLI_PARSER_H
#define CLI_PARSER_H

#include <iostream>
#include <map>
#include <sstream>
#include <string>

class cli_parser
{
public:
  cli_parser(int argc, char *argv[])
  {
    auto i = argc - 1;
    while (i)
    {
      std::string key, val;
      if (argv[i][0] == '-' && !isdigit(argv[i][1])) // flaga
      {
        key = argv[i];
        val = "1";
      }
      else // zmienna
      {
        key = argv[i - 1];
        val = argv[i--];
      }
      --i;
      m[key] = val;
    }
  }

  template <typename T = int>
  T get(std::string key)
  {
    if (!m.count("-" + key))
    {
      throw std::runtime_error("missing key: -" + key);
    }
    std::stringstream ss(m["-" + key]);
    T result;
    ss >> result;
    return result;
  }

  template <typename T>
  T get(std::string key, T default_value)
  {
    return m.count("-" + key) ? get<T>(key) : default_value;
  }

  // -fconcepts ?
  // decltype(auto) get(std::string key, auto default_value)
  // {
  //   if (m.count("-" + key))
  //     return get<decltype(default_value)>(key);
  //   else
  //     return default_value;
  // }

  bool get_flag(std::string key) { return get<bool>(key, false); }

private:
  std::map<std::string, std::string> m;
};

#endif // CLI_PARSER_H