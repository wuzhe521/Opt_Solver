# MATRIX CLASS

Matrix 是开源项目GONS的核心类，用于矩阵运算。模块具有轻量化，效率高和使用简单的特点。  
为保持兼容性和用户友好性，Matrix类的设计参考了MATLAB的矩阵操作方式。  
在实际的使用中，Matrix类可以直接使用，也可以对Matrix类进行封装，以便于使用。  
Matrix 类的实现不依赖于任何第三方库。
Matrix 类以头文件的方式实现，无需额外为Matrix类处理编译问题，方便用户IOTO - "Include One Time Only" 的使用。  
[[Basic Matrix Operation refered from MATLAB](https://ww2.mathworks.cn/help/matlab/math/matrices-in-the-matlab-environment.html)]  


![Linux build status](https://github.com/autodiff/autodiff/workflows/linux/badge.svg?branch=master)
![Windows build status](https://github.com/autodiff/autodiff/workflows/windows/badge.svg?branch=master)
![Coverage](https://codecov.io/gh/autodiff/autodiff/branch/master/graph/badge.svg)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![GitHub release](https://img.shields.io/github/release/autodiff/autodiff.svg)](https://github.com/autodiff/autodiff/releases)
[![GitHub issues](https://img.shields.io/github/issues/autodiff/autodiff.svg)](https://github.com/autodiff/autodiff/issues)


#矩阵类
========================================================================================
## 矩阵的创建
1. 默认创建
```cpp
Matrix<float, 3, 3> mat;
```
2. 拷贝
```cpp
Matrix<float, 3, 3> mat_copy(mat);
```
## 矩阵值的赋值和获取
1. 按索引赋值
```cpp
Matrix<float, 3, 3> mat;
mat(0, 0) = 1;
mat(1, 1) = 2;
mat(2, 2) = 3;
```
2. 初始化赋值
```cpp
Matrix<float, 3, 3> mat{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
```
3. 获取
```cpp
float value = mat(0, 0);
```
## 矩阵的行列数
1. 获取行数
```cpp
int row = mat.Row();
int col = mat.Col();
```
## 矩阵的加法和减法
矩阵的加减法是逐个元素执行的，或者说是按元素执行的
```cpp
Matrix<float, 3, 3> A{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
Matrix<float, 3, 3> B{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
Matrix<float, 3, 3> C = A + B;
C.Print("A + B = \n");
```
>A + B =  
>  2 2 2  
>  2 2 2  
>  2 2 2  

```cpp
Matrix<float, 3, 3> A{{2, 2, 2}, {2, 2, 2}, {2, 2, 2}};
Matrix<float, 3, 3> B{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
Matrix<float, 3, 3> C = A - B;
C.Print("A - B = \n");
```
>A - B =  
>  1 1 1  
>  1 1 1  
>  1 1 1  


## 矩阵的乘法
如果 A 的列维度等于 B 的行维度，或者其中一个矩阵为标量，则可定义矩阵乘积 C = AB。如果 A 为 m×p 且 B 为 p×n，则二者的乘积 C 为 m×n。矩形矩阵乘法必须满足维度兼容性条件,如果维度不一致，会产生错误。乘法不能以相反的顺序执行。
```cpp
gons::Matrix<float, 2, 3> M1{{1, 2, 3}, {4, 5, 6}};
gons::Matrix<float, 3, 2> M2{{7, 8}, {9, 10}, {11, 12}};
gons::Matrix<float, 2, 2> M3 = M1 * M2;
M3.Print("M1 * M2 = \n");
```

>M1 * M2 =   
>58 64   
>139 154 

```cpp
gons::Matrix<float, 2, 3> M1{{1, 2, 3}, {4, 5, 6}};
gons::Matrix<float, 4, 4> M4 = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
auto M5 = M1 * M4; // This should trigger a CHECK error
```
>ERROR: Matrix multiplication dimension mismatch: 3 != 4

## 矩阵的转置
对于实矩阵，转置运算是对矩阵A中的元素 aij 和 aji 进行交换。数学中常用A' 或 A^T 来表示
```cpp
gons::Matrix<float, 2, 3> M{{1, 2, 3}, {4, 5, 6}};
M1.Print("M : \n");
gons::Matrix<float, 3, 2> MT = M.Transpose();
MT.Print("M' : \n");
```
>M :   
>1 2 3   
>4 5 6   
>M' :   
>1 4   
>2 5   
>3 6

## 单位矩阵
单位矩阵是对角线上元素都是1，其余元素都是0的方阵， 假设单位矩阵是Ann， n是矩阵的维度，A(i, i) = 1, A(i, j) = 0 i nequal j 。 单位矩阵的生成右两种方式，通过已有方阵eye() 生成，和通过特殊函数 Identity()生成。
```cpp
gons::Matrix<float, 4, 4> M4 = {
      {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
auto M4_I = M4.eye();
M4_I.Print("M4_I eye() : \n");
```
> M4_I eye() :  
> 1 0 0 0   
> 0 1 0 0   
> 0 0 1 0   
> 0 0 0 1  
```cpp
auto M5_I = gons::Matrix<float, 4, 4>::Identify();
M5_I.Print("M5_I Identify() :");
```
> M5_I Identify() :  
> 1 0 0 0   
> 0 1 0 0   
> 0 0 1 0   
> 0 0 0 1  
## 矩阵的逆
矩阵的逆是矩阵的逆矩阵，也称为伴生矩阵。
如果矩阵 A 为非奇异方阵（非零行列式），则方程 AX = I 和 XA = I 具有相同的解 X。此解称为 A 的逆矩阵，表示为 $A^{-1}$ 。  
```cpp
gons::Matrix<float, 2, 2> M4 = {
      {1, 2}, {3, 4}};
auto M4_I = M4.Inverse();
M4_I.Print("M4_I Inverse() : \n");
```
> M4_I Inverse() :  
> -2 1   
> 1.5 -0.5  

如果矩阵 A 为奇异矩阵（零行列式），则 A 的逆不存在。Matrux 类中，如果矩阵 A 为奇异矩阵，则 Inverse() 函数会返回一个系统错误。

```cpp
gons::Matrix<float, 3, 3> M3 = {
      {1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
auto M3_I = M3.Inverse();
M3_I.Print("M3_I Inverse() : \n");
```
>ERROR: Matrix is singular.

## 向量
向量是矩阵的特殊情况，向量是具有一行或一列的矩阵。  

### 向量的创建
gons也支持向量的创建。在gons中，向量是矩阵的特殊情况，向量是具有一行或一列的矩阵。默认是列向量。
```cpp
gons::Vector<float, 3> vec ={1, 2, 3};  
```
### 向量的加法和减法
向量的加减法是逐个元素执行的，或者说是按元素执行的
```cpp
gons::Vector<float, 3> vec2 ={4, 5, 6};  
auto vec3 = vec + vec2;
vec3.Print("vec3 = vec + vec2 : \n");
```
> vec3 = vec + vec2 :  
> 5 7 9  
## 向量值的获取和设置
向量值的获取和矩阵稍有不同，向量值的获取和设置是通过[]索引来实现的。
```cpp
gons::Vector<float, 3> vec ={1, 2, 3};  
vec.Print("vec : \n");
vec(0) = 10;
vec.Print("vec(0) = 10 : \n");
```
> vec :  
> 1 2 3  
> vec(0) = 10 :  
> 10 2 3  
