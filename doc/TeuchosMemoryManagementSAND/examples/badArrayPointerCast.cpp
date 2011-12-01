#include <iostream>


class A {
  char dummy_;
  char *char_ptr_;
public:
  A() : dummy_('a'), char_ptr_(&dummy_) {}
  void incrementA() { ++(*char_ptr_); }
  char get_dummy() const { return dummy_; }
};


class B : public A {
  int size_;
public:
  B() : size_(0) {}
};



void goodFoo(const int n) {
  std::cout << "\ngoodFoo("<<n<<"):\n";
  A *a_array = new A[n];
  for (int i = 0; i < n; ++i) {
    a_array[i].incrementA(); // KABOMMMMM!
    std::cout << "a_array["<<i<<"] = '"<<a_array[i].get_dummy()<<"'\n";
  }
  delete [] a_array;
}


void badFoo(const int n) {
  std::cout << "\nbadFoo("<<n<<"):\n";
  B *b_array = new B[n];
  A *a_array = b_array; // Compiles just fine :-(
  for (int i = 0; i < n; ++i) {
    a_array[i].incrementA(); // KABOMMMMM!
    std::cout << "a_array["<<i<<"] = '"<<a_array[i].get_dummy()<<"'\n";
  }
  delete [] b_array;
}


int main(){
  goodFoo(2);
  badFoo(2);
  return 0;
}
