# Inertial_navigation_Quadcopter
四轴飞行器中惯性导航算法C语言实现

利用加速度计陀螺仪地磁计采集的数据，用四元数的方法计算出欧拉角，其中偏航角用互补滤波计算

`IMUupdate`:由九轴输入计算出欧拉角

`Hcm_correction1` :地磁计矫正

移植输入九轴传感器的数字可算出欧拉角
