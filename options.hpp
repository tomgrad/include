#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <fstream>
#include <map>
#include <sstream>
#include <string>

class options
{
public:
  options(int argc, char *argv[])
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
      m[key.substr(1)] = val;
    }
  }

  void parsefile(const std::string &filename, const bool overwrite_cmdline = false)
  {
    std::ifstream f(filename);

    while (f.good())
    {
      std::string line;
      std::getline(f, line);
      if (line[0] == '#')
        continue;
      auto pos = line.find('=');
      const auto key = line.substr(0, pos);
      const auto val = line.substr(++pos);
      if (!m.count(key) || overwrite_cmdline)
        m[key] = val;
    }
    f.close();
  }

  template <typename T = int>
  T get(std::string key)
  {
    if (!m.count(key))
    {
      throw std::runtime_error("missing key: -" + key);
    }
    std::stringstream ss(m[key]);
    T result;
    ss >> result;
    return result;
  }

  template <typename T>
  T get(std::string key, T default_value)
  {
    return m.count(key) ? get<T>(key) : default_value;
  }

  bool get_flag(std::string key) { return get<bool>(key, false); }

private:
  std::map<std::string, std::string> m;
};

#endif // OPTIONS_HPP