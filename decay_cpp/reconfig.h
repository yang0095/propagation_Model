#ifndef RECONFIG_H
#define RECONFIG_H
#include <iostream>
#include <string>
#include<map>
#include<vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <cstdio>  // 用于 printf
#ifdef _WIN32
#include <windows.h> // 用于 GetCurrentDirectory（Windows）
#else
#include <unistd.h> // 用于 getcwd（Linux）
#include <dirent.h>  // POSIX 目录操作
#endif
using namespace  std;

class  ReConfig{
public:
      std::map<std::string, std::map<std::string, std::string> >settings_;

    string getTime(int hour);

    string getDirPath(string path_str);

    // 检测系统字节序
    bool isLittleEndian();

    // 将4字节数据从大端转换为本地字节序
    double toFloatBigEndian(const char* data);


    // 读取大端存储的二进制文件
    std::vector<double> readEra15Hm(const std::string& filePath, int month,int ilat, int ilong );

    void Trim(std::string & str);


    bool IsSpace(char c);


    bool AnalyseLine(const std::string & line, std::string& section, std::string & key, std::string & value);



    bool ReadConfig(const std::string & filename);


    std::string ReadString(const char* section, const char* item, const char* default_value);


};

#endif // RECONFIG_H
