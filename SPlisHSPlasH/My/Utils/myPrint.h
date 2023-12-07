
#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip> 

template<typename T>
void printVectorNohup(T contents, int maxTimes=999, std::string fileName="", int precision=8, std::string msg="")
{
    static int times = 0;
    times++;
    if (times>maxTimes)
        return;
    
    if(fileName=="")
    {
        std::cout<<msg;
        for(auto x:contents)
            std::cout<<std::fixed <<std::setprecision(precision)<<x<<"\t";
        std::cout<<"\n";
        return;
    }
    else
    {
        static std::ofstream fout;
        fout.open(fileName, std::ios::app);
        fout<<msg;
        for(const auto& x:contents)
            fout<<std::fixed <<std::setprecision(precision)<<x<<"\t";
        fout<<"\n";
        fout.close();
        return;
    }
}

template<typename T>
void printVec(T contents, int precision = 8)
{
    for(auto x:contents)
        std::cout<<std::fixed <<std::setprecision(precision)<<x<<"\t";
    std::cout<<"\n";
    return;
}

template<typename T>
void printScalarNohup(T contents, int maxTimes = 100, std::string fileName="debugOutput.txt",int precision=16)
{
    static int times = 0;
    times++;
    if (times>maxTimes)
        return;
    
    static std::ofstream fout;
    fout.open(fileName, std::ios::app);
    fout<<std::fixed <<std::setprecision(precision)<<contents<<"\n";
    fout.close();
    return;
}

/**
 * @brief 保存向量场到txt
 * 
 * @tparam T 
 * @param fileName 文件名 
 * @param content 要打印的场
 * @param precision 精度（默认小数点后8位数）
 */
template<typename T>
void printVectorField(std::string fileName, T content,  size_t precision=8)
{
    std::ofstream f;
    f.open(fileName);
    for(const auto& x:content)
    {
        for(const auto& xx:x)
            f<<std::fixed <<std::setprecision(precision)<<xx<<"\t";
        f<<"\n";
    } 
    f.close();
}

/**
 * @brief 保存标场到txt
 * 
 * @tparam T 
 * @param fileName 文件名 
 * @param content 要打印的场
 * @param precision 精度（默认小数点后8位数）
 */
template<typename T>
void printScalarField(std::string fileName, T content,  size_t precision=8)
{
    std::ofstream f;
    f.open(fileName);
    for(const auto& x:content)
    {
        f<<std::fixed <<std::setprecision(precision)<<x<<"\n";
    } 
    f.close();
}

#define echo(content) {std::cout<<(#content)<<": "<<content<<std::endl;}
#define echoVec(content) {std::cout<<(#content)<<": ";  printVec(content);}