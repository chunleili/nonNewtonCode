#pragma once
#include "SPlisHSPlasH/Common.h"

//用户交互（键盘鼠标）的中介类，用于传递和处理数据
struct Interactive
{
    Vector3r mouse_pos; //鼠标位置
    enum KEY{left=0, right, up, down, forward, backward};
    Vector3r* p_rb_pos = nullptr; //刚体质心位置的指针

    //a singleton method to get the object
    static Interactive& get_inter()
    {
        static Interactive inter;
        return inter;
    }
    
    //把mouse_pos从外界传递给Interactive内部
    void get_mouse_pos(const Vector3r& rhs)
    {
        mouse_pos[0] = rhs[0];
        mouse_pos[1] = rhs[1];
        mouse_pos[2] = rhs[2];
        printf("Mouse world space pos: (%.3f, %.3f, %.3f)\n", mouse_pos[0],mouse_pos[1],mouse_pos[2]);
    }

    //获取键盘的输入：从GUI\OpenGL\MiniGL.cpp MiniGL::keyboard
    void get_key_input(enum KEY input)
    {
        const float move_speed = 0.1;

        if(input == KEY::down)
        {
            if(p_rb_pos != nullptr)
            {
                (*p_rb_pos)[1] -= move_speed;
                printf("down!\n");
            }
        }
        else if(input == KEY::up)
        {
            if(p_rb_pos != nullptr)
            {
                (*p_rb_pos)[1] += move_speed;
                printf("up!\n");
            }
        }
        else if(input == KEY::left)
        {
            if(p_rb_pos != nullptr)
            {
                (*p_rb_pos)[0] -= move_speed;
                printf("left!\n");
            }
        }
        else if (input == KEY::right)
        {
            if (p_rb_pos != nullptr)
            {
                (*p_rb_pos)[0] += move_speed;
                printf("right!\n");
            }
        }
        else if (input == KEY::forward)
        {
            if (p_rb_pos != nullptr)
            {
                (*p_rb_pos)[2] -= move_speed;
                printf("forward!\n");
            }
        }
        else if (input == KEY::backward)
        {
            if (p_rb_pos != nullptr)
            {
                (*p_rb_pos)[2] += move_speed;
                printf("backward!\n");
            }
        }
    }

    //获取刚体的控制权。在PBDWrapper.cpp中调用
    void set_rb_pos(Vector3r& rb_pos)
    {
        // 获取并设定位置为鼠标点击位置 FIXME:
        // (rb_pos) = mouse_pos; 

        //用wasdfb来移动
        // 获取刚体质心指针
        p_rb_pos = &rb_pos; 
        // std::cout<< "rb_pos: "<< (rb_pos)<<"\n";
    }
};