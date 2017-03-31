#pragma once

#include <set>
#include <vector>
#include <list>
#include <map>
#include <cstdint>

namespace IKD{
    template <class T>
    class cell;

    template<class T>
    class tree_object{
    public:
        cell<T> *cell = nullptr;
        T *object = nullptr;
        tree_object<T> *prev = nullptr;
        tree_object<T> *next = nullptr;
        const int id;

        tree_object(int id = 0) : id(id){}
        ~tree_object() = default;

        bool remove(){
            if(!cell){
                return false;
            }
            if(!cell->on_remove(this)){
                return false;
            }

            if(prev){
                prev->next = next;
            }
            if(next){
                next->prev = prev;
            }
            prev = nullptr;
            next = nullptr;
            cell = nullptr;
            return true;
        }

    };

    template < class T >
    class collision_list {
    public:
        static const std::size_t realloc_size = 1000;

        collision_list(){
            vec.reserve(realloc_size);
        }

        ~collision_list() = default;

        size_t size() {
            return vec.size();
        }

        T **root() {
            return &vec.front();
        }

        void reset() {
            vec.resize(0);
        }

        void write(T *obj1, T *obj2){
            if(vec.size() == vec.capacity()){
                vec.reserve(vec.capacity() + realloc_size);
            }
            vec.push_back(obj1);
            vec.push_back(obj2);
        }

    private:
        std::vector<T*> vec;
    };

    #define CLINER4TREEMANAGER_MAXLEVEL        9
    template <class T>
    class liner_for_tree_manager{
    private:
        cell<T> **cell_array;
        unsigned int pow[CLINER4TREEMANAGER_MAXLEVEL + 1];
        double width;
        double height;
        double left;
        double top;
        double unit_width;
        double unit_height;
        std::uint32_t cell_num;
        unsigned int level;
        collision_list<T> collision_list_instance;

    public:
        liner_for_tree_manager(){
            level = 0;
            width = 0.0;
            height = 0.0;
            left = 0.0;
            top = 0.0;
            unit_width = 0.0;
            unit_height = 0.0;
            cell_num = 0;
            cell_array = nullptr;

            pow[0] = 1;
            for(int i = 1; i < CLINER4TREEMANAGER_MAXLEVEL + 1; i++){
                pow[i] = pow[i - 1] * 4;
            }
        }

        ~liner_for_tree_manager(){
            for(std::uint32_t i = 0; i < cell_num; i++){
                if(cell_array[i] != nullptr){
                    delete cell_array[i];
                }
            }
            delete[] cell_array;
        }

        bool init(unsigned int Level, double left, double top, double right, double bottom){
            if(Level >= CLINER4TREEMANAGER_MAXLEVEL){
                return false;
            }

            cell_num = (pow[Level + 1] - 1) / 3;
            cell_array = new cell<T>*[cell_num]{ nullptr };

            left = left;
            top = top;
            width = right - left;
            height = bottom - top;
            unit_width = width / (1 << Level);
            unit_height = height / (1 << Level);

            level = Level;

            return true;
        }

        bool register_object(double left, double top, double right, double bottom, tree_object<T> *oft){
            std::uint32_t elem = get_morton_number(left, top, right, bottom);
            if(elem < cell_num){
                if(!cell_array[elem]){
                    create_new_cell(elem);
                }
                return cell_array[elem]->push(oft);
            }
            return false;
        }

        std::uint32_t get_all_collision_list(collision_list<T> *&list){
            collision_list_instance.reset();
            if(cell_array[0] == nullptr){
                return 0;
            }

            std::list<T*> stack;
            get_collision_list(0, stack);
            list = &collision_list_instance;
            return (std::uint32_t)collision_list_instance.size();
        }



    private:
        bool get_collision_list(std::uint32_t elem, std::list<T*> &stack){
            tree_object<T> *oft1 = cell_array[elem]->get_first_object();
            while(oft1 != 0){
                tree_object<T> *oft2 = oft1->next;
                while(oft2 != 0){
                    collision_list_instance.write(oft1->object, oft2->object);
                    oft2 = oft2->next;
                }
                for(std::list<T*>::iterator it = stack.begin(); it != stack.end(); it++){
                    collision_list_instance.write(oft1->object, *it);
                }
                oft1 = oft1->next;
            }

            bool child_flag = false;
            std::uint32_t obj_num = 0;
            std::uint32_t next_elem;
            for(std::uint32_t i = 0; i < 4; i++){
                next_elem = elem * 4 + 1 + i;
                if(next_elem < cell_num && cell_array[elem * 4 + 1 + i]){
                    if(!child_flag){
                        oft1 = cell_array[elem]->get_first_object();
                        while(oft1 != 0){
                            stack.push_back(oft1->object);
                            obj_num++;
                            oft1 = oft1->next;
                        }
                    }
                    child_flag = true;
                    get_collision_list(elem * 4 + 1 + i, stack);
                }
            }

            if(child_flag){
                for(std::uint32_t i = 0; i < obj_num; i++){
                    stack.pop_back();
                }
            }

            return true;
        }


        bool create_new_cell(std::uint32_t elem){
            while(!cell_array[elem]){
                cell_array[elem] = new cell<T>;
                elem = (elem - 1) >> 2;
                if(elem >= cell_num){
                    break;
                }
            }
            return true;
        }

        std::uint32_t get_morton_number(double left, double top, double right, double bottom){
            std::uint32_t lt = get_point_element(left, top);
            std::uint32_t rb = get_point_element(right, bottom);

            std::uint32_t def = rb ^ lt;
            unsigned int hi = 0;
            for(unsigned int i = 0; i < level; i++)
            {
                std::uint32_t Check = (def >> (i * 2)) & 0x3;
                if(Check != 0){
                    hi = i + 1;
                }
            }
            std::uint32_t space_num = rb >> (hi * 2);
            std::uint32_t add_num = (pow[level - hi] - 1) / 3;
            space_num += add_num;

            if(space_num > cell_num){
                return 0xffffffff;
            }

            return space_num;
        }

        std::uint32_t bit_separate32(std::uint32_t n){
            n = (n | (n << 8)) & 0x00ff00ff;
            n = (n | (n << 4)) & 0x0f0f0f0f;
            n = (n | (n << 2)) & 0x33333333;
            return (n | (n << 1)) & 0x55555555;
        }

        std::uint16_t get_2d_morton_number(std::uint16_t x, std::uint16_t y){
            return static_cast<std::uint16_t>(bit_separate32(x) | (bit_separate32(y) << 1));
        }

        std::uint32_t get_point_element(double pos_x, double pos_y){
            return get_2d_morton_number(
                static_cast<std::uint16_t>((pos_x - left) / unit_width),
                static_cast<std::uint16_t>((pos_y - top) / unit_height)
            );
        }
    };

    template <class T>
    class cell{
    private:
        tree_object<T> *lastest_oft = nullptr;

    public:
        cell() = default;

        ~cell(){
            if(lastest_oft){
                reset_link(lastest_oft);
            }
        }

        void reset_link(tree_object<T> *oft){
            if(oft->next) {
                reset_link(oft->next);
            }
            oft = nullptr;
        }

        bool push(tree_object<T> *oft){
            if(oft == nullptr){
                return false;
            }
            if(oft->cell == this){
                return false;
            }
            if(lastest_oft == nullptr){
                lastest_oft = oft;
            }else{
                oft->next = lastest_oft;
                lastest_oft->prev = oft;
                lastest_oft = oft;
            }
            oft->cell = this;
            return true;
        }

        tree_object<T> *get_first_object(){
            return lastest_oft;
        }

        bool on_remove(tree_object<T> *remove_obj)
        {
            if(lastest_oft == remove_obj){
                if(lastest_oft != nullptr){
                    lastest_oft = lastest_oft->next;
                }
            }
            return true;
        }
    };
}