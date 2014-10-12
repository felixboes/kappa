// The following example demonstrates the execution of a private member function unsing templates.
// The example is due to
//
//     http://bloglitb.blogspot.de/2010/07/access-to-private-members-thats-easy.html
//
#include<iostream>

template<typename Tag>
struct result {
    typedef typename Tag::type type;
    static type ptr;
};

template<typename Tag>
typename result<Tag>::type result<Tag>::ptr;

template<typename Tag, typename Tag::type p>
struct rob : result<Tag> {
    struct filler {
        filler()
        {
            result<Tag>::ptr = p;
        }
    };
    static filler filler_obj;
};

template<typename Tag, typename Tag::type p>
typename rob<Tag, p>::filler rob<Tag, p>::filler_obj;

class A {
private:
    void f() {
        std::cout << "proof!" << std::endl;
    }
};

struct Af {
    typedef void(A::*type)();
};
template class rob<Af, &A::f>;

int main() {
  A a;
  (a.*result<Af>::ptr)();
  return 0;
}
