#include <vector>
#include <iostream>

template<class T>
class BlockTransformerBase {
public:
  virtual ~BlockTransformerBase() {};
  virtual void transform(const std::vector<T> &a, std::vector<T> &b) const =0;
};

template<class T>
class AddIntoTransformer : public BlockTransformerBase <T> {
public:
  virtual void transform(const std::vector<T> &a, std::vector<T> &b) const {
    for (int i = 0; i < a.size(); ++i)
      b[i] += a[i];
  }
};

int main(){
  const int xx[]={1,2,3,4,5};
  std::vector <int> x(xx,xx+5);
  const int yy[]={6,7,8,9,10};
  std::vector <int> y(yy,yy+5);
  
  AddIntoTransformer <int>().transform(x,y);
  
  for(std::vector <int>::iterator z=y.begin();z!=y.end();z++){
    std::cout << *z << std::endl;
  }
  return 0;
}
