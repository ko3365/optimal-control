# Control Optimization of the System with Disturbance and Sensor Noise using LQR and Kalman Filter

### Given the input force _u_ applied to the cart, the disturbance _d_ acting on the pendulum tip, and sensor noise _v_, design an optimal controller that estimates the state and minimizes the cost so that the closed-loop system converges and the inverted pendulum tip is stabilized.

<p align="center">
  <img 
    width="250"
    height="275"
    src="images/system_pic.PNG"
  >
</p>

The physical system is modeled with three state variables: pendulum angle, angular velocity, cart velocity. The system is linearlized assuming small pendulum angle.
M = 2kg, m = 1kg, l = 4m, g = 9.8m/s^2

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}&space;x&space;=&space;\begin{bmatrix}&space;\theta\\\dot{\theta}\\\dot{x_c}&space;\end{bmatrix},&space;\dot{x}&space;=&space;Ax&plus;Bu&plus;Ld,&space;A&space;=&space;\begin{bmatrix}&space;0&space;&&space;1&space;&&space;0&space;\\&space;\frac{(M&plus;m)g}{Ml}&space;&&space;0&space;&&space;0&space;\\&space;\frac{-mg}{M}&0&0\end{bmatrix},&space;B&space;=&space;\begin{bmatrix}&space;0\\\frac{-1}{Ml}\\\frac{1}{M}\end{bmatrix},&space;L&space;=&space;\begin{bmatrix}&space;0&space;\\&space;\frac{1}{ml}\\0\end{bmatrix}&space;}">

## Step 1. Design Feedback Controller (LQR)
Let's say that we have a sensor for the cart velocity and we want to have our system optimally converge to the static equilibrium. The cost function is defined to be:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}&space;J&space;=\int_{0}^{\infty}&space;\dot{x_c}^2&space;&plus;&space;(u-u_o)^2&space;\&space;\&space;dt}">

Define the perturbation variables as:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}&space;\xi&space;=&space;x-x_o,&space;\&space;\&space;\mu&space;=&space;u-u_o}">

Then, the Q and R matrices can be calculated:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}&space;C&space;=&space;\begin{bmatrix}&space;0&&space;0&space;&1\end{bmatrix},&space;\&space;\&space;\&space;\dot{x_c}^2&space;=&space;\xi^TC^TC\xi,&space;\&space;\&space;\&space;&space;(u-u_o)^2&space;=&space;\mu^TR\mu,&space;\&space;\&space;\&space;Q&space;=&space;C^TC,&space;\&space;\&space;R&space;=&space;1}">

Three necessary and sufficient conditions for the existence of an optimal control are checked, and then lqr function is used to obtain the controller gain.

1. (A,B) stabilizable -> PBH Rank Test: rank(ctrb(A,B)) = 3
2. Q and R non-singular positive matrix
3. (Q,A) has no unobservable mode on imaginary axis: rank(obsv(C,A)) = 3

```Matlab
[Ko,~,~] = lqr(A,B,Q,R);
Ko = -Ko;
```

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}&space;\mu&space;=&space;K_o\xi&space;\rightarrow&space;u&space;=&space;K_ox&plus;K_1d}">

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}&space;K_o&space;=&space;\begin{bmatrix}&space;69.4076&space;&&space;37.5637&space;&&space;1.0000\end{bmatrix},&space;K_1&space;=&space;\begin{bmatrix}K_o&space;&-I\end{bmatrix}\begin{bmatrix}A&B\\C&0\end{bmatrix}^{-1}&space;\begin{bmatrix}L\\0\end{bmatrix}&space;=&space;6.0824}">

## Step 2. Design State Estimator (Kalman Filter)

Define the augmented system and introduce disturbance _d_ and sensor noise _v_
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}&space;\dot{\chi}&space;=&space;&space;\mathcal{A}\chi&plus;&space;\mathbb{B}u&space;&plus;&space;\mathcal{B}f,&space;\&space;\&space;\chi&space;=&space;\begin{bmatrix}&space;x\\d&space;\end{bmatrix},&space;\&space;\&space;\mathcal{A}=&space;\begin{bmatrix}&space;A&space;&&space;L&space;\\&space;0&space;&&space;0&space;&space;\end{bmatrix},&space;\&space;\&space;&space;\mathbb{B}&space;=&space;\begin{bmatrix}&space;B\\0&space;\end{bmatrix},&space;\&space;\&space;\mathcal{B}&space;=&space;\begin{bmatrix}&space;0\\I&space;\end{bmatrix},&space;\&space;\&space;&space;f&space;=&space;d_o\delta(t)}">

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{\dot{\hat{\chi}}&space;=&space;&space;\mathcal{A}\hat{\chi}&plus;&space;\mathbb{B}u&plus;F(\mathcal{C}\hat{\chi}-y),&space;\&space;\&space;\&space;y=&space;\mathcal{C}\chi&plus;v}}">

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{\dot{e}&space;=&space;\dot{\chi}-\dot{\hat{\chi}}&space;=&space;\mathcal{A}e&plus;\mathcal{B}f&plus;\underbrace{F(\underbrace{\mathcal{C}e&plus;v)}_\text{$\eta$}}_\text{$u$},&space;\&space;\&space;\&space;\&space;z&space;=&space;e}}">

Generalized plant for H2 optimal Control:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{\begin{bmatrix}&space;\dot{e}\\z\\&space;\eta&space;\end{bmatrix}&space;=&space;\begin{bmatrix}&space;&space;&space;\mathcal{A}&space;&&space;\mathcal{B}&space;&&space;0&space;&&space;I&space;\\&space;&space;&space;I&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;&space;&space;\mathcal{C}&space;&&space;0&space;&&space;I&&space;0&space;\end{bmatrix}&space;\begin{bmatrix}&space;e\\f\\v\\u\end{bmatrix}&space;}">

Applying the Duality Property of LQR and Kalman Filter, we can obtain the Kalman Gain F:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{W&space;=&space;\mathcal{B}\mathcal{B}^T&space;\geq&space;0,&space;\&space;\&space;V&space;=&space;I&space;>0&space;}">

```Matlab
[F,~,~] = lqr(Ak',Ck',W,V);
F = -F';
```

## Step 3. Simulation (Basic Model)

With the controller gain _K_ and the Kalman gain _F_, we can now define the closed loop system and simulate the basic model (sensor noise ignored)

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{u&space;=&space;K\hat{\chi},&space;\&space;\&space;\dot{\chi}&space;=&space;&space;\mathcal{A}\chi&plus;&space;\mathbb{B}u&space;&plus;&space;\mathcal{B}f,&space;\&space;\&space;\dot{\hat{\chi}}&space;=&space;&space;\mathcal{A}\hat{\chi}&plus;&space;\mathbb{B}u&plus;F(\mathcal{C}\hat{\chi}-y)&space;}">

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{\begin{bmatrix}\dot{\chi}\\\dot{\hat{\chi}}\\y\end{bmatrix}&space;=&space;\begin{bmatrix}&space;&space;&space;\mathcal{A}&space;&&space;\mathbb{B}K&space;&&space;\mathcal{B}&space;&&space;0&space;\\&space;&space;&space;-F\mathcal{C}&space;&&space;\mathcal{A}&plus;\mathbb{B}K&plus;F\mathcal{C}&space;&&space;0&space;&&space;-F&space;\\&space;&space;&space;\mathcal{C}&space;&&space;0&space;&&space;0&&space;I\end{bmatrix}&space;\begin{bmatrix}\chi\\\hat{\chi}\\f\\v\end{bmatrix}&space;}">

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{A_{cl}&space;=&space;\begin{bmatrix}\mathcal{A}&space;&&space;\mathbb{B}K\\-F\mathcal{C}&space;&&space;\mathcal{A}&plus;\mathbb{B}K&plus;F\mathcal{C}\end{bmatrix},&space;\&space;\&space;\chi(0)&space;=&space;\begin{bmatrix}0\\0\\0\\1\end{bmatrix},&space;\&space;\&space;&space;\hat{\chi}(0)=\begin{bmatrix}0\\0\\0\\0\end{bmatrix}&space;}">

Since there is no input to the basic model (no sensor noise), we can use init function in Matlab to simulate the system:

```Matlab
model_basic = ss(Acl,[],[],[]);
[~,t,x] = init(model_basic,x0,20);
```

<p align="center">
  <img 
    width="700"
    height="500"
    src="images/basic_design.PNG"
  >
</p>

## Step 4. Simulation (Additional Tuning Parameter)
Introduce tuning parameter so that our system can be robust with the presence of the sensor noise:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{u&space;\rightarrow&space;\alpha&space;u,&space;\&space;\&space;\&space;v&space;\rightarrow&space;v/\beta&space;}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{R=I&space;\rightarrow&space;\alpha^2,&space;\&space;\&space;V&space;=&space;I&space;\rightarrow&space;\beta^2}">

With new parameters, controller gain and Kalman gain need to be re-obtained

```Matlab
[Ko,~,~] = lqr(A,B,Q,alpha^2);
Ko = -Ko;
[F,~,~] = lqr(Ak',Ck',W,beta^2);
F = -F';
```
Tuning Method: 

-Increasing alpha will slower the state convergence and decreasing alpha will fasten the state convergence.

-Increasing beta will reduce the noise sensitivity. Tune alpha to meet the convergence requirements and tune beta to meet the noise sensitivity requirements.

System is also defined so that it contains the tuning parameters:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{\begin{bmatrix}\dot{x}\\\dot{\hat{\chi}}\\z1&space;\triangleq&space;\dot{x_c}\\z2&space;\triangleq&space;\alpha&space;u\end{bmatrix}&space;=&space;\begin{bmatrix}&space;&space;&space;A&space;&&space;BK&space;&&space;L&space;&&space;0&space;\\&space;&space;&space;-FC&space;&&space;\mathcal{A}&plus;\mathbb{B}K&plus;F\mathcal{C}&space;&&space;0&space;&&space;-\beta&space;F&space;\\&space;&space;&space;C&space;&&space;0&space;&&space;0&&space;0\\&space;&space;&space;0&space;&&space;\alpha&space;K&space;&&space;0&space;&&space;0\end{bmatrix}\begin{bmatrix}x\\&space;\chi\\\begin{pmatrix}d\\v/\beta\end{pmatrix}\end{bmatrix}}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{A_{cl}&space;=&space;\begin{bmatrix}A&space;&&space;BK\\-FC&space;&&space;\mathcal{A}&plus;\mathbb{B}K&plus;F\mathcal{C}\end{bmatrix},&space;\&space;\&space;B_{cl}&space;=&space;\begin{bmatrix}&space;L&space;&&space;0\\&space;0&space;&&space;-\beta&space;F\end{bmatrix},&space;\&space;\&space;C_{cl}&space;=&space;\begin{bmatrix}C&space;&&space;0\\&space;0&space;&&space;\alpha&space;K\end{bmatrix},&space;\&space;\&space;D_{cl}&space;=&space;\begin{bmatrix}0&&space;0\\0&0\end{bmatrix}}">

Assuming that the sensor noise has a sinusoidal property:
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}{x_c(0)=0,\&space;\&space;x(0)=0,\&space;\&space;\hat{\chi}(0)=0,\&space;\&space;d(t)=1,&space;\&space;\&space;v(t)/\beta=0.01sin(20\pi&space;t)/\beta}">

```Matlab
model_robust = ss(Acl,Bcl,Ccl,Dcl);
t = linspace(0,20,1000);
u1 = ones(1,length(t));
u2 = 0.01*sin(20*pi*t)/beta;
[y,t,x] = lsim(model_robust,[u1' u2'],t,x0);
```
<p align="center">
  <img 
    width="900"
    height="400"
    src="images/design_tuning.PNG"
  >
</p>

Due to the sensor noise, there is a fluctuation in controller input; however, the fluctuation is minimized with the tuning parameters. The figure on the right shows the zoom-ed in graph of input _u_, and we can see that the amplitude after 14 seconds (max-min) are less than 0.4. The following design specifications are satisfied:

(i) The cart-pendulum converges to static equilibrium.

(ii) The cart-pendulum comes to a stop with the minimum displacement. (LQR)

(iii) The system is robust against the sensor noise. (Input amplitude after 14sec is less than 0.04)


