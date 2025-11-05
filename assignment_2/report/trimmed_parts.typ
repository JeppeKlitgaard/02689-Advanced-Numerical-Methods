= Tymek 1.b

Consider
$
  1/r partial / (partial r) (r (partial phi)/ (partial r) (r, theta)) + 1/r^2 (partial^2 phi) / (partial theta ^ 2) (r, theta) = 0 \
  r (partial phi) / (partial r) (r, theta) + r^2 (partial^2 phi)/ (partial r^2) (r, theta) + (partial^2 phi) / (partial theta ^ 2) (r, theta) = 0 \
  phi(r_1, theta) = g_1(theta); quad phi(r_2, theta) = g_2(theta) \
  phi(r, 0) = phi(r, 2 pi) = g_3(r)
$
Exact solution
$
  phi(r, theta) = V_infinity (r + r_1^2 / r) cos(theta)
$
Consider
$
  r = (r_2 - r_1)/2 z_r + (r_2 + r_1)/2, quad partial / (partial r) = 2/(r_2 - r_1) partial / (partial z_r), quad partial^2 / (partial r^2) = 4/(r_2 - r_1)^2 partial / (partial z_r)
$
Define
$
  alpha_1 = (r_2 - r_1)/2, alpha_2 = (r_2 + r_1)/2, quad r = alpha_1 z_r + alpha_2, quad partial / (partial r) = 1/alpha_1 partial / (partial z_r), quad partial^2 / (partial r^2) = 1/alpha_1^2 partial / (partial z_r)
$
Consider
$
  theta = (2 pi)/2 z_theta + (2 pi)/2 = pi z_theta + pi, quad partial / (partial theta) = 1/pi partial / (partial z_theta), quad partial^2 / (partial theta^2) = 1/pi^2 partial^2 / (partial z_theta^2)
$
Consider
$
    r (partial phi) / (partial r) + r^2 (partial^2 phi)/ (partial r^2) + (partial^2 phi) / (partial theta ^ 2) = 0 \

    (alpha_1 z_r + alpha_2) 1/alpha_1 partial / (partial z_r) phi + (alpha_1 z_r + alpha_2)^2 1/alpha_1^2 partial^2 / (partial z_r^2) phi + 1/pi^2 partial^2 / (partial z_theta^2) phi = 0 \

      (z_r + alpha_2 / alpha_1) partial / (partial z_r) phi + (z_r + alpha_2 / alpha_1)^2 partial^2 / (partial z_r^2) phi + 1/pi^2 partial^2 / (partial z_theta^2) phi = 0
$
Solution becomes
$
  phi(z_r, z_theta) = (alpha_1 z_r + alpha_2 + r_1^2 / (alpha_1 z_r + alpha_2)) cos(pi z_theta + pi)
$
Collocation method discretization
$
  D_x = "diag"(beta_1, ..., beta_N) D_r + "diag"(beta_1^2, ..., beta_N^2) D_r^2\
  D_y = 1/pi^2 D_theta^2\
  L_N = "kron"(I, D_x) + "kron"(D_y, I)
$
We define the new coordinate system $(z_r, z_theta) = (x, y)$, we use the row-wise ordering.

Two dimensional Vandermonde matrices
$
  U = V_x hat(U) V_y^T
$
