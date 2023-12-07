#include <string>

#include "Partio.h"
#include "PartioSingleton.h"

using namespace Partio;


// 测试一下在随便一个函数和随便一个类中使用PartioSingleton::getCurrent()是否能够正常工作
void func()
{
    auto *d = PartioSingleton::getCurrent();
    std::cout << "-----------I am in func" << std::endl;

    d->print();
}

class someClass
{
public:
    someClass()
    {
        auto *d = PartioSingleton::getCurrent();
        std::cout << "----------I am in someClass" << std::endl;
        d->print();
    }
};

int main()
{
    auto *d = PartioSingleton::getCurrent();
    // d->read("ice-cream.bhclassic");
    // d->print();

    d->read("test.bgeo");
    // d->print();

    func();

    someClass a;

    return 0;
}