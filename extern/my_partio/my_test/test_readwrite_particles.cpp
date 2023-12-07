#include <vector>
#include <string>
#include <array>
#include "Partio.h"
#include "WriteToPartio.h"
using namespace Partio;
int main()
{
    //makeData
    float spacing = 0.1;
    std::vector<std::array<float,3>> a(100);
    for (size_t i = 0; i < 100; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            a[i][j] = i * spacing ;
        }
    }
    
    // call writePartio in PartioIO.h
    const std::string outputFile("E:/Dev/SPlisHSPlasH/extern/my_partio/my_test/test_out.bgeo.gz");
    writePartio(outputFile, a);
    std::cout<<"write to "<<outputFile<<"\n";
}