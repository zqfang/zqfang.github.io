# Derivative, Gradient, Jacobian, Hessian, Laplacian


Just some basic notations

## Derivative

$$
f^{\prime} (x) = \frac {df(x)} {dx}
$$

## Gradient 

Generalize the derivative to the multivariate functions.   
The first order derivative of a multivariate functions.

$$
\nabla f=\left[\frac{\partial f\left(x_{1}, x_{2}, x_{3}\right)}{\partial x_{1}}, \frac{\partial f\left(x_{1}, x_{2}, x_{3}\right)}{\partial x_{2}}, \frac{\partial f\left(x_{1}, x_{2}, x_{3}\right)}{\partial x_{3}}\right]
$$

## Jacobian

a generalization of the derivate operator to the vector-valued functions

$$
J=\left(\begin{array}{cccc}
\frac{\partial f_{1}}{\partial x_{1}} & \frac{\partial f_{1}}{\partial x_{2}} & \cdots & \frac{\partial f_{1}}{\partial x_{n}} \cr
\frac{\partial f_{2}}{\partial x_{1}} & \frac{\partial f_{2}}{\partial x_{2}} & \cdots & \frac{\partial f_{2}}{\partial x_{n}} \cr
\vdots & \vdots & \ddots & \vdots \cr
\frac{\partial f_{m}}{\partial x_{1}} & \frac{\partial f_{m}}{\partial x_{2}} & \cdots & \frac{\partial f_{m}}{\partial x_{n}}
\end{array}\right)
$$


## Hessian

the second order derivative of a multivariate function

$$
H=\left(\begin{array}{cccc}
\frac{\partial^{2} f}{\partial x_{1}^{2}} & \frac{\partial^{2} f}{\partial x_{1} \partial x_{2}} & \cdots & \frac{\partial^{2} f}{\partial x_{1} \partial x_{n}} \cr
\frac{\partial^{2} f}{\partial x_{2} \partial x_{1}} & \frac{\partial^{2} f}{\partial x_{2}^{2}} & \cdots & \frac{\partial^{2} f}{\partial x_{2} \partial x_{n}} \cr
\vdots & \vdots & \ddots & \vdots \cr
\frac{\partial^{2} f}{\partial x_{n} \partial x_{1}} & \frac{\partial^{2} f}{\partial x_{n} \partial x_{2}} & \cdots & \frac{\partial^{2} f}{\partial x_{n}^{2}}
\end{array}\right)
$$

## Laplacian
The trace of the Hessian matrix is known as the Laplacian operator denoted by $\nabla^2$:

$$
\nabla^2 f = \operatorname{trace}(H) = \frac{\partial^{2} f}{\partial x_{1}^{2}}+\frac{\partial^{2} f}{\partial x_{2}^{2}}+\cdots+\frac{\partial^{2} f}{\partial x_{n}^{2}}
$$

## Laplace's equation
The second order partial differential equation

In rectangular coordinates, 

$$
\nabla^{2} f=\frac{\partial^{2} f}{\partial x^{2}}+\frac{\partial^{2} f}{\partial y^{2}}+\frac{\partial^{2} f}{\partial z^{2}}=0
$$

In cylindrical coordinates,

$$
\nabla^{2} f=\frac{1}{r} \frac{\partial}{\partial r}\left(r \frac{\partial f}{\partial r}\right)+\frac{1}{r^{2}} \frac{\partial^{2} f}{\partial \phi^{2}}+\frac{\partial^{2} f}{\partial z^{2}}=0
$$

In spherical coordinates, using the $(r, \theta, \varphi)$ convention, 

$$
\nabla^{2} f=\frac{1}{r^{2}} \frac{\partial}{\partial r}\left(r^{2} \frac{\partial f}{\partial r}\right)+\frac{1}{r^{2} \sin \theta} \frac{\partial}{\partial \theta}\left(\sin \theta \frac{\partial f}{\partial \theta}\right)+\frac{1}{r^{2} \sin ^{2} \theta} \frac{\partial^{2} f}{\partial \varphi^{2}}=0
$$

More generally, in curvilinear coordinates,

$$
\nabla^{2} f=\frac{\partial}{\partial \xi^{j}}\left(\frac{\partial f}{\partial \xi^{k}} g^{k j}\right)+\frac{\partial f}{\partial \xi^{j}} g^{j m} \Gamma_{m n}^{n}=0
$$

or


$$
\nabla^{2} f=\frac{1}{\sqrt{|g|}} \frac{\partial}{\partial \xi^{i}}\left(\sqrt{|g|} g^{i j} \frac{\partial f}{\partial \xi^{j}}\right)=0, \quad\left(g=\operatorname{det}\left\{g_{i j}\right\}\right)
$$


[Laplace's equation](https://en.wikipedia.org/wiki/Laplace%27s_equation)
