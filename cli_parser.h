#ifndef CLI_PARSER_H
#define CLI_PARSER_H

#include <iostream>
#include <map>
#include <sstream>
#include <string>

class cli_parser {
public:
  cli_parser(int argc, char *argv[]) {
    auto i = argc - 1;
    while (i) {
      std::string key, val;
      if (argv[i][0] == '-') // flaga
      {
        key = argv[i];
        val = "1";
      } else // zmienna
      {
        key = argv[i - 1];
        val = argv[i--];
      }
      --i;
      m[key] = val;
    }
  }

  template <typename T> T get(std::string key) {
    if (!m.count("-" + key)) {
      std::cerr << "wrong key: " << key << std::endl;
      exit(1);
    }
    std::stringstream ss(m["-" + key]);
    T result;
    ss >> result;
    return result;
  }

  template <typename T> T get(std::string key, T default_value) {
    if (m.count("-" + key))
      return get<T>(key);
    else
      return default_value;
  }

  bool get_flag(std::string key) { return get<bool>(key, false); }

  // size_t is_key(std::string key) { return m.count("-" + key); }

private:
  std::map<std::string, std::string> m;
};

#endif // CLI_PARSER_H