#include 
...

struct mesh{
        args...;
};

typedef vector<vector<struct mesh>> grid;

struct init_cond{
        thermostats;
        pressiostats;
};

class eqdiff
{
        public:
                void launch_expe();
                void iteration();
                void remesh();
                void eqdiff();
                attr gets();
                void sets;

        private:
                unsigned long int nb_of_ite;
                unsigned int, size_x, size_y;
                grid datas;
                struct init_cond init_conditions;
};      



