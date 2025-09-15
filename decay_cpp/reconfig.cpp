#include"reconfig.h"

string ReConfig::getTime(int hour){
    int minval=0,diff=24;
    for(int k=0;k<=18;k+=6){
        if(abs(k-hour)<diff){
            minval=k;
        }
    }
    string time=to_string(minval);

    if(minval<10){
        time="0"+time;
    }

    return time;
}


std::string ReConfig::getDirPath(string path_str) {
    char buffer[1024]; // 创建一个缓冲区来存储路径

    // 根据操作系统选择合适的函数
#ifdef _WIN32
        // 获取当前工作目录（Windows）
    DWORD length = GetCurrentDirectoryA(sizeof(buffer), buffer);
    if (length == 0) {
        std::cerr << "获取当前目录失败" << std::endl;
        return ""; // 返回空字符串表示错误
    }
#else
        // 获取当前工作目录（Linux 和其他 POSIX 系统）
    if (getcwd(buffer, sizeof(buffer)) == nullptr) {
        std::cerr << "获取当前目录失败" << std::endl;
        return ""; // 返回空字符串表示错误
    }
#endif

    std::string configPath = buffer; // 获取当前路径
    //printf("Cfg=%s\n", configPath.c_str()); // 打印当前路径

    // 假设 ReadConfig() 和 ReadString() 已经定义
    // 读取配置文件
//    bool ret = ReadConfig(configPath + "/path.ini");
//    if (ret == false) {
//        printf("读取配置文件出错, Cfg=%s\n", "path.ini"); // 打印错误信息
//        return ""; // 返回空字符串表示错误
//    }

    // 从配置文件中读取指定的路径
    const char*pathStr=path_str.c_str();
    std::string dirPath = ReadString("path", pathStr, "");
    //std::cout << dirPath << std::endl; // 打印读取到的路径

    return dirPath; // 返回读取到的路径
}


// 检测系统字节序
bool ReConfig::isLittleEndian() {
    uint16_t num = 1;
    return *(reinterpret_cast<uint8_t*>(&num)) == 1;
}

// 将4字节数据从大端转换为本地字节序
double ReConfig::toFloatBigEndian(const char* data) {
    uint32_t value;
    std::memcpy(&value, data, sizeof(value)); // Copy 4 bytes into uint32_t

    // 如果本地字节序是小端，则需要进行字节序转换
    if (isLittleEndian()) {
        value = ((value & 0xFF000000) >> 24) |
                ((value & 0x00FF0000) >> 8) |
                ((value & 0x0000FF00) << 8) |
                ((value & 0x000000FF) << 24);
    }

    float result;
    std::memcpy(&result, &value, sizeof(result));
    return (double)result;
}


// 读取大端存储的二进制文件
std::vector<double> ReConfig::readEra15Hm(const std::string& filePath, int month,int ilat, int ilong ) {
    std::ifstream file(filePath, std::ios::binary);
    if (!file) {
        std::cerr << "无法打开文件: " << filePath << std::endl;
        return {};
    }

    std::vector<double> data(22);
    std::vector<char> buffer(4); // 4 bytes for float
    for(int i=1;i<=22;++i){
        for(int j=1;j<=12;++j){
            file.read(buffer.data(), 4);
            if (!file) {
                std::cerr << "读取文件时出错。" << std::endl;
                return {};
            }
            if(j==month){
                data[i-1] = toFloatBigEndian(buffer.data());
            }
        }
    }
    return data;
}


void ReConfig::Trim(std::string & str)
{
    if (str.empty())
    {
        return;
    }
    int i, start_pos, end_pos;
    for (i = 0; i < str.size(); ++i) {
        if (!IsSpace(str[i])) {
            break;
        }
    }
    if (i == str.size())
    {
        str = "";
        return;
    }
    start_pos = i;
    for (i = str.size() - 1; i >= 0; --i) {
        if (!IsSpace(str[i])) {
            break;
        }
    }
    end_pos = i;
    str = str.substr(start_pos, end_pos - start_pos + 1);
}

bool ReConfig::IsSpace(char c)
{
    if (' ' == c || '\t' == c)
        return true;
    return false;
}


bool ReConfig::AnalyseLine(const std::string & line, std::string& section, std::string & key, std::string & value)
{
    if (line.empty())
        return false;
    int start_pos = 0, end_pos = line.size() - 1, pos, s_startpos, s_endpos;
    if ((pos = line.find("#")) != -1)
    {
        if (0 == pos)
        {
            return false;
        }
        end_pos = pos - 1;
    }
    if (((s_startpos = line.find("[")) != -1) && ((s_endpos = line.find("]"))) != -1)
    {
        section = line.substr(s_startpos + 1, s_endpos - 1);
        return true;
    }
    std::string new_line = line.substr(start_pos, start_pos + 1 - end_pos);
    if ((pos = new_line.find('=')) == -1)
        return false;
    key = new_line.substr(0, pos);
    value = new_line.substr(pos + 1, end_pos + 1 - (pos + 1));
    Trim(key);
    if (key.empty()) {
        return false;
    }
    Trim(value);
    if ((pos = value.find("\r")) > 0)
    {
        value.replace(pos, 1, "");
    }
    if ((pos = value.find("\n")) > 0)
    {
        value.replace(pos, 1, "");
    }
    return true;
}


bool ReConfig::ReadConfig(const std::string & filename)
{
    settings_.clear();
    std::ifstream infile(filename.c_str());//构造默认调用open,所以可以不调用open

    if (!infile) {
        return false;
    }
    std::string line, key, value, section;
    std::map<std::string, std::string> k_v;
    std::map<std::string, std::map<std::string, std::string> >::iterator it;
    while (getline(infile, line))
    {
        if (AnalyseLine(line, section, key, value))
        {
            it = settings_.find(section);
            if (it != settings_.end())
            {
                k_v[key] = value;
                it->second = k_v;
            }
            else
            {
                k_v.clear();
                settings_.insert(std::make_pair(section, k_v));
            }
        }
        key.clear();
        value.clear();
    }
    infile.close();
    return true;
}

std::string ReConfig::ReadString(const char* section, const char* item, const char* default_value)
{
    std::string tmp_s(section);
    std::string tmp_i(item);
    std::string def(default_value);
    std::map<std::string, std::string> k_v;
    std::map<std::string, std::string>::iterator it_item;
    std::map<std::string, std::map<std::string, std::string> >::iterator it;
    it = settings_.find(tmp_s);
    if (it == settings_.end())
    {

        return def;
    }
    k_v = it->second;
    it_item = k_v.find(tmp_i);
    if (it_item == k_v.end())
    {

        return def;
    }
    return it_item->second;
}
