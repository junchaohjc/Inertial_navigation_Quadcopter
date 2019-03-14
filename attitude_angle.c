/****************************************************************************
* 名    称：attitude_angle

* 说    明：计算欧拉角
            利用加速度计陀螺仪地磁计采集的数据
						用四元数的方法计算出欧拉角
						其中偏航角用互补滤波计算
* 作		者 :junchaohjc
****************************************************************************/ 

#include "attitude_angle.h"
#include "MPU6050.h"
#include "math.h"
#include "HMC5883L.h"
#include "usart.h"

IMU_DATA GyroFinal;             //读取处理后化为rad/s单位的角速度量
IMU_DATA AccFinal;             //读取处理后化为m/s^2单位的加速度量
IMU_DATA Geomagnetic;           //地磁计
S_FLOAT_ANGLE Q_ANGLE;          //姿态欧拉角
S_FLOAT_ANGLE Q_ANGLE_OFFSET;          //姿态欧拉角偏移 //利用此值改变航向角期望，使每次飞行开始航向角期望皆为0 //


float gravity;

static float q0 = 1.0, q1 = 0.0, q2 = 0.0, q3 = 0.0;
// int16_t hcm_x, hcm_y, hcm_z;   				//读取地磁计初值，椭圆矫正时会用，故用全局变量

float ROTATE_XY[4]={ 1.0170,0.0016,0.0016,1.0001};//{0.0018,0,0,0.0018};
float ROTATE_YZ[4]={ 0.9339,-0.0682,-0.0682,0.9296};//{0.0020,-0.0001,-0.0001,0.0020};
/****************************************************************************
* 名    称：IMUupdate
* 功    能：由九轴输入计算出欧拉角
* 入口参数：加速度，角速度，地磁场
* 出口参数：
* 说    明：内部使用了静态变量
* 调用方法：无 
****************************************************************************/ 
void IMUupdate(float gx, float gy, float gz, float ax, float ay, float az)
{   
//  static float exInt = 0.0f, eyInt = 0.0f, ezInt = 0.0f;
	
// float delta_2=0;
// float delta_4=0;
  float norm=1.0;
	float gravityerror=1.0;
  //float hx, hy, hz, bx, bz;
  float M11,M21,M31,M13,M23,M33;  //计算出来的参考系Z轴（重力）在机体坐标系中的向量
	float NX_X,NX_Y,NX_Z;
  float GZ_ex=0, GZ_ey=0, GZ_ez=0,MX_ex=0,MX_ey=0,MX_ez=0;//参考系与计算出的载体系之间的误差

  float q0q0 = q0*q0;
  float q0q1 = q0*q1;
  float q0q2 = q0*q2;
  float q0q3 = q0*q3;
  float q1q1 = q1*q1;
  float q1q2 = q1*q2;
  float q1q3 = q1*q3;
  float q2q2 = q2*q2;
  float q2q3 = q2*q3;
  float q3q3 = q3*q3;
	
  norm = sqrt(ax*ax + ay*ay +az*az);       // 测量正常化 把加速度计的三维向量转成单位向量。
// 	printf(" x %f y %f z %f\n",ax,ay,az);

// 	gravityerror=fabs(gravity-norm);
// 	printf(" x %f y %f z %f\n",norm ,gravity,gravityerror);
// 	if(gravityerror>1) {gravityerror=1;}
// 	gravityerror=1-gravityerror*gravityerror; //参数两用，误差越大，重力矫正时误差乘的系数越小
  ax = ax /norm;
  ay = ay /norm;
  az = az /norm;
//     	printf(" x %f y %f z %f\n",norm ,gravity,gravityerror);
  M13= 2*(q1q3 - q0q2);						// 计算参考坐标系Z轴在机体坐标系下的向量					
  M23= 2*(q0q1 + q2q3);
  M33= q0q0-q1q1-q2q2+q3q3;
	
	M11=q0q0+q1q1-q2q2-q3q3;         //计算参考坐标系x轴在机体坐标系下的向量
	M21=2*(q1q2-q0q3);
	M31=2*(q1q3+q0q2);
	
// 	//首先根据地磁向量与重力向量计算出一个在参考坐标系中代表x轴的向量M*G
// 	NX_X=Geomagnetic.Y*M33-Geomagnetic.Z*M23;
// 	NX_Y=Geomagnetic.Z*M13-Geomagnetic.X*M33;
// 	NX_Z=Geomagnetic.X*M23-Geomagnetic.Y*M13;
// 	//单位化
// 	norm=sqrt(NX_X*NX_X+NX_Y*NX_Y+NX_Z*NX_Z);
// 	NX_X=NX_X/norm;
// 	NX_Y=NX_Y/norm;
// 	NX_Z=NX_Z/norm;
// // 做计算出的在机体坐标系中的x轴向量乘实际值，注意顺序，
// 	MX_ex=M31*NX_Y-M21*NX_Z;
// 	MX_ey=M11*NX_Z-M31*NX_X;
// 	MX_ez=M21*NX_X-M11*NX_Y;

// 	MX_ex=0;
// 	MX_ey=0;
// 	MX_ez=0;

//这是把四元数换算成《方向余弦矩阵》中的第三列的三个元素。
//根据余弦矩阵和欧拉角的定义，地理坐标系的重力向量，转到机体坐标系，正好是这三个元素。
//所以这里的vx\y\z，其实就是当前的欧拉角（即四元数）的机体坐标参照系上，换算出来的重力单位向量。


  // Error is sum of cross product between estimated direction and measured direction of field vectors
  GZ_ex = (ay*M33 - az*M23);                       //注意向量叉乘顺序    					
  GZ_ey = (az*M13 - ax*M33);
  GZ_ez = (ax*M23 - ay*M13);
	
	
//axyz是机体坐标参照系上，加速度计测出来的重力向量，也就是实际测出来的重力向量。
//axyz是测量得到的重力向量，vxyz是陀螺积分后的姿态来推算出的重力向量，它们都是机体坐标参照系上的重力向量。
//那它们之间的误差向量，就是陀螺积分后的姿态和加计测出来的姿态之间的误差。
//向量间的误差，可以用向量叉积（也叫向量外积、叉乘）来表示，exyz就是两个重力向量的叉积。
//这个叉积向量仍旧是位于机体坐标系上的，而陀螺积分误差也是在机体坐标系，而且叉积的大小与陀螺积分误差成正比，正好拿来纠正陀螺。（你可以自己拿东西想象一下）由于陀螺是对机体直接积分，所以对陀螺的纠正量会直接体现在对机体坐标系的纠正。

/*  
  exInt = exInt + ex * Ki ;					//计算和应用积分反馈			 
  eyInt = eyInt + ey * Ki ;
  ezInt = ezInt + ez * Ki ;


  gx = gx + Kp*ex + exInt;					//校正陀螺仪测量值	   用叉积误差来做PI修正陀螺零偏							
  gy = gy + Kp*ey + eyInt;
  gz = gz + Kp*ez + ezInt;		  */ 

//  halfT=GET_NOWTIME();			  				//两次计算的时间间隔，单位秒

//   gx = gx + gravityerror*GZ_ex+0.2*MX_ex; 					//校正陀螺仪测量值	   用叉积误差来做PI修正陀螺零偏							
//   gy = gy + gravityerror*GZ_ey+0.2*MX_ey; 
//   gz = gz + gravityerror*GZ_ez+0.2*MX_ez;	 

// 	
  gx = gx + GZ_ex; 					//校正陀螺仪测量值	   用叉积误差来做PI修正陀螺零偏							
  gy = gy + GZ_ey; 
  gz = gz + GZ_ez;
/* 
  delta_2=(2*halfT*gx)*(2*halfT*gx)+(2*halfT*gy)*(2*halfT*gy)+(2*halfT*gz)*(2*halfT*gz);	
  delta_4=delta_2*delta_2;

  q0 = (1-delta_2/8+delta_4/384)*q0 + (-q1*gx - q2*gy - q3*gz)*halfT*(0.5-delta_2/48);			// 整合四元数率	 四元数微分方程	四元数更新算法，四阶毕卡法
  q1 = (1-delta_2/8+delta_4/384)*q1 + (q0*gx + q2*gz - q3*gy)*halfT*(0.5-delta_2/48);
  q2 = (1-delta_2/8+delta_4/384)*q2 + (q0*gy - q1*gz + q3*gx)*halfT*(0.5-delta_2/48);
  q3 = (1-delta_2/8+delta_4/384)*q3 + (q0*gz + q1*gy - q2*gx)*halfT*(0.5-delta_2/48);	*/		 

 
//   delta_2=(2*halfT*gx)*(2*halfT*gx)+(2*halfT*gy)*(2*halfT*gy)+(2*halfT*gz)*(2*halfT*gz);	
//  
//   q0 = (1-delta_2/8)*q0 + (-q1*gx - q2*gy - q3*gz)*halfT;			// 整合四元数率	 四元数微分方程	四元数更新算法，二阶毕卡法
//   q1 = (1-delta_2/8)*q1 + (q0*gx + q2*gz - q3*gy)*halfT;
//   q2 = (1-delta_2/8)*q2 + (q0*gy - q1*gz + q3*gx)*halfT;
//   q3 = (1-delta_2/8)*q3 + (q0*gz + q1*gy - q2*gx)*halfT;			 
 												   
  q0 = q0 + (-q1*gx - q2*gy - q3*gz)*halfT;							// 整合四元数率	 四元数微分方程	四元数更新算法，一阶龙库法
  q1 = q1 + (q0*gx + q2*gz - q3*gy)*halfT;
  q2 = q2 + (q0*gy - q1*gz + q3*gx)*halfT;         
  q3 = q3 + (q0*gz + q1*gy - q2*gx)*halfT;			   	   
//   mmx= q0 + (-q1*gx - q2*gy - q3*gz)*halfT;							// 整合四元数率	 四元数微分方程	四元数更新算法，一阶龙库法
//   mmy= q1 + (q0*gx + q2*gz - q3*gy)*halfT;
//   mmz= q2 + (q0*gy - q1*gz + q3*gx)*halfT;         
//   mmf= q3 + (q0*gz + q1*gy - q2*gx)*halfT;		

  norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);		// 正常化四元
  q0 = q0 / norm;
  q1 = q1 / norm;
  q2 = q2 / norm;
  q3 = q3 / norm;
// printf(" %f; %f; %f; %f\n",q0,q1,q2,q3);
//   Q_ANGLE.Pitch = asin(-2*q1*q3 + 2*q0*q2)*57.3; 					//pitch
//   Q_ANGLE.Roll = atan2(2*q2*q3 + 2*q0*q1, -2*q1*q1 - 2*q2*q2+1)*57.3; 	// roll
//   Q_ANGLE.Yaw = -atan2(2*q1*q2 + 2*q0*q3, -2*q2*q2 - 2*q3*q3+1)*57.3; // yaw
//   //转换为欧拉角,两种计算结果均有误差，周期中间量不可以准确计算
  Q_ANGLE.Roll =-asin(M13)*57.3; // pitch
	Q_ANGLE.Pitch = atan2(M23,M33)* 57.3; 	// roll
  Q_ANGLE.Yaw = -atan2(2*q1q2 + 2*q0q3,M11)* 57.3; // yaw
//   		printf("Q_ANGLE.Roll %f  Q_ANGLE.Pitch %f  AngEYaw % f\n",Q_ANGLE.Roll,Q_ANGLE.Pitch,Q_ANGLE.Yaw);

	
}




/****************************************************************************
* 名    称：Hcm_correction1()
* 功    能：地磁计矫正
* 入口参数：无
* 出口参数：无
* 说    明：
* 调用方法：无 
****************************************************************************/ 


void  Hcm_correction()
{
	u8  hbuffer[8];
 	int16_t hcm_x, hcm_y, hcm_z;   				//读取地磁计初值
	
  IICreadBytes(HMC58X3_ADDR,HMC58X3_R_XM,6,hbuffer);
	hcm_x=(((int16_t)hbuffer[0]) << 8) | hbuffer[1];
	hcm_z=(((int16_t)hbuffer[2]) << 8) | hbuffer[3];
  hcm_y=(((int16_t)hbuffer[4]) << 8) | hbuffer[5];                      //手动根据MATLAB椭圆矫正算法来矫正地磁计的周围磁场强度
	

	hcm_x= (hcm_x- HMC_OFFSET_X); 
 	hcm_y= (hcm_y- HMC_OFFSET_Y);
	hcm_z= (hcm_z- HMC_OFFSET_Z); //减去一个偏移值，然后乘以一个缩放量，以Y轴为基准 
	
	Geomagnetic.X = ROTATE_XY[0]*hcm_x+ROTATE_XY[1]*hcm_y; 
 	Geomagnetic.Y = ROTATE_XY[2]*hcm_x+ROTATE_XY[3]*hcm_y;
	Geomagnetic.Z = ROTATE_YZ[2]*hcm_y+ROTATE_YZ[3]*hcm_z;
}
