#include "pch.h"
#include "tensor.h"

int main(int argc, char** argv)
{
    using namespace std;
    size_t n1{4}, n2{4};
    size_t n1m{n1-1}, n2m{n2-1};
    size_t n1p{n1+1}, n2p{n2+1};
    Tensor<double> 
        e1(n1p), 
        e2(n2p);
    double e1_start{0.}, e1_end{1.};
    double e2_start{0.}, e2_end{1.};
    double de1 {(e1_end-e1_start)/static_cast<double>(n1)};
    double de2 {(e1_end-e1_start)/static_cast<double>(n2)};

    for (size_t i = 0; i < e1.size(); ++i)
        e1[i] = e1_start + static_cast<double>(i)*de1;

    for (size_t j = 0; j < e2.size(); ++j)
        e2[j] = e2_start + static_cast<double>(j)*de2;

#ifdef DEBUG
    cout << "----------grid-----------" << endl;
    cout << e1 << endl;
    cout << e2 << endl;    
#endif

    Tensor<Tensor<double>> 
        u(n1p, Tensor<double>(n2p))     , v(n1p, Tensor<double>(n2p)),
        u0(n1p, Tensor<double>(n2p))    , v0(n1p, Tensor<double>(n2p)),
        he1(n1p, Tensor<double>(n2p))   , he2(n1p, Tensor<double>(n2p)),
        phi(n1p, Tensor<double>(n2p))   , p(n1p, Tensor<double>(n2p)),
        psi(n1p, Tensor<double>(n2p))   , vor(n1p, Tensor<double>(n2p)),
        stream(n1p, Tensor<double>(n2p)), vort(n1p, Tensor<double>(n2p));
        
    u = 0.;
                       


    return 0;
}


