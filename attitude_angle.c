/****************************************************************************
* ��    �ƣ�attitude_angle

* ˵    ��������ŷ����
            ���ü��ٶȼ������ǵشżƲɼ�������
						����Ԫ���ķ��������ŷ����
						����ƫ�����û����˲�����
* ��		�� :junchaohjc
****************************************************************************/ 

#include "attitude_angle.h"
#include "MPU6050.h"
#include "math.h"
#include "HMC5883L.h"
#include "usart.h"

IMU_DATA GyroFinal;             //��ȡ�����Ϊrad/s��λ�Ľ��ٶ���
IMU_DATA AccFinal;             //��ȡ�����Ϊm/s^2��λ�ļ��ٶ���
IMU_DATA Geomagnetic;           //�شż�
S_FLOAT_ANGLE Q_ANGLE;          //��̬ŷ����
S_FLOAT_ANGLE Q_ANGLE_OFFSET;          //��̬ŷ����ƫ�� //���ô�ֵ�ı亽���������ʹÿ�η��п�ʼ�����������Ϊ0 //


float gravity;

static float q0 = 1.0, q1 = 0.0, q2 = 0.0, q3 = 0.0;
// int16_t hcm_x, hcm_y, hcm_z;   				//��ȡ�شżƳ�ֵ����Բ����ʱ���ã�����ȫ�ֱ���

float ROTATE_XY[4]={ 1.0170,0.0016,0.0016,1.0001};//{0.0018,0,0,0.0018};
float ROTATE_YZ[4]={ 0.9339,-0.0682,-0.0682,0.9296};//{0.0020,-0.0001,-0.0001,0.0020};
/****************************************************************************
* ��    �ƣ�IMUupdate
* ��    �ܣ��ɾ�����������ŷ����
* ��ڲ��������ٶȣ����ٶȣ��شų�
* ���ڲ�����
* ˵    �����ڲ�ʹ���˾�̬����
* ���÷������� 
****************************************************************************/ 
void IMUupdate(float gx, float gy, float gz, float ax, float ay, float az)
{   
//  static float exInt = 0.0f, eyInt = 0.0f, ezInt = 0.0f;
	
// float delta_2=0;
// float delta_4=0;
  float norm=1.0;
	float gravityerror=1.0;
  //float hx, hy, hz, bx, bz;
  float M11,M21,M31,M13,M23,M33;  //��������Ĳο�ϵZ�ᣨ�������ڻ�������ϵ�е�����
	float NX_X,NX_Y,NX_Z;
  float GZ_ex=0, GZ_ey=0, GZ_ez=0,MX_ex=0,MX_ey=0,MX_ez=0;//�ο�ϵ������������ϵ֮������

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
	
  norm = sqrt(ax*ax + ay*ay +az*az);       // ���������� �Ѽ��ٶȼƵ���ά����ת�ɵ�λ������
// 	printf(" x %f y %f z %f\n",ax,ay,az);

// 	gravityerror=fabs(gravity-norm);
// 	printf(" x %f y %f z %f\n",norm ,gravity,gravityerror);
// 	if(gravityerror>1) {gravityerror=1;}
// 	gravityerror=1-gravityerror*gravityerror; //�������ã����Խ����������ʱ���˵�ϵ��ԽС
  ax = ax /norm;
  ay = ay /norm;
  az = az /norm;
//     	printf(" x %f y %f z %f\n",norm ,gravity,gravityerror);
  M13= 2*(q1q3 - q0q2);						// ����ο�����ϵZ���ڻ�������ϵ�µ�����					
  M23= 2*(q0q1 + q2q3);
  M33= q0q0-q1q1-q2q2+q3q3;
	
	M11=q0q0+q1q1-q2q2-q3q3;         //����ο�����ϵx���ڻ�������ϵ�µ�����
	M21=2*(q1q2-q0q3);
	M31=2*(q1q3+q0q2);
	
// 	//���ȸ��ݵش��������������������һ���ڲο�����ϵ�д���x�������M*G
// 	NX_X=Geomagnetic.Y*M33-Geomagnetic.Z*M23;
// 	NX_Y=Geomagnetic.Z*M13-Geomagnetic.X*M33;
// 	NX_Z=Geomagnetic.X*M23-Geomagnetic.Y*M13;
// 	//��λ��
// 	norm=sqrt(NX_X*NX_X+NX_Y*NX_Y+NX_Z*NX_Z);
// 	NX_X=NX_X/norm;
// 	NX_Y=NX_Y/norm;
// 	NX_Z=NX_Z/norm;
// // ����������ڻ�������ϵ�е�x��������ʵ��ֵ��ע��˳��
// 	MX_ex=M31*NX_Y-M21*NX_Z;
// 	MX_ey=M11*NX_Z-M31*NX_X;
// 	MX_ez=M21*NX_X-M11*NX_Y;

// 	MX_ex=0;
// 	MX_ey=0;
// 	MX_ez=0;

//���ǰ���Ԫ������ɡ��������Ҿ����еĵ����е�����Ԫ�ء�
//�������Ҿ����ŷ���ǵĶ��壬��������ϵ������������ת����������ϵ��������������Ԫ�ء�
//���������vx\y\z����ʵ���ǵ�ǰ��ŷ���ǣ�����Ԫ�����Ļ����������ϵ�ϣ����������������λ������


  // Error is sum of cross product between estimated direction and measured direction of field vectors
  GZ_ex = (ay*M33 - az*M23);                       //ע���������˳��    					
  GZ_ey = (az*M13 - ax*M33);
  GZ_ez = (ax*M23 - ay*M13);
	
	
//axyz�ǻ����������ϵ�ϣ����ٶȼƲ����������������Ҳ����ʵ�ʲ����������������
//axyz�ǲ����õ�������������vxyz�����ݻ��ֺ����̬����������������������Ƕ��ǻ����������ϵ�ϵ�����������
//������֮�������������������ݻ��ֺ����̬�ͼӼƲ��������̬֮�����
//������������������������Ҳ�������������ˣ�����ʾ��exyz�����������������Ĳ����
//�����������Ծ���λ�ڻ�������ϵ�ϵģ������ݻ������Ҳ���ڻ�������ϵ�����Ҳ���Ĵ�С�����ݻ����������ȣ����������������ݡ���������Լ��ö�������һ�£����������ǶԻ���ֱ�ӻ��֣����Զ����ݵľ�������ֱ�������ڶԻ�������ϵ�ľ�����

/*  
  exInt = exInt + ex * Ki ;					//�����Ӧ�û��ַ���			 
  eyInt = eyInt + ey * Ki ;
  ezInt = ezInt + ez * Ki ;


  gx = gx + Kp*ex + exInt;					//У�������ǲ���ֵ	   �ò���������PI����������ƫ							
  gy = gy + Kp*ey + eyInt;
  gz = gz + Kp*ez + ezInt;		  */ 

//  halfT=GET_NOWTIME();			  				//���μ����ʱ��������λ��

//   gx = gx + gravityerror*GZ_ex+0.2*MX_ex; 					//У�������ǲ���ֵ	   �ò���������PI����������ƫ							
//   gy = gy + gravityerror*GZ_ey+0.2*MX_ey; 
//   gz = gz + gravityerror*GZ_ez+0.2*MX_ez;	 

// 	
  gx = gx + GZ_ex; 					//У�������ǲ���ֵ	   �ò���������PI����������ƫ							
  gy = gy + GZ_ey; 
  gz = gz + GZ_ez;
/* 
  delta_2=(2*halfT*gx)*(2*halfT*gx)+(2*halfT*gy)*(2*halfT*gy)+(2*halfT*gz)*(2*halfT*gz);	
  delta_4=delta_2*delta_2;

  q0 = (1-delta_2/8+delta_4/384)*q0 + (-q1*gx - q2*gy - q3*gz)*halfT*(0.5-delta_2/48);			// ������Ԫ����	 ��Ԫ��΢�ַ���	��Ԫ�������㷨���ĽױϿ���
  q1 = (1-delta_2/8+delta_4/384)*q1 + (q0*gx + q2*gz - q3*gy)*halfT*(0.5-delta_2/48);
  q2 = (1-delta_2/8+delta_4/384)*q2 + (q0*gy - q1*gz + q3*gx)*halfT*(0.5-delta_2/48);
  q3 = (1-delta_2/8+delta_4/384)*q3 + (q0*gz + q1*gy - q2*gx)*halfT*(0.5-delta_2/48);	*/		 

 
//   delta_2=(2*halfT*gx)*(2*halfT*gx)+(2*halfT*gy)*(2*halfT*gy)+(2*halfT*gz)*(2*halfT*gz);	
//  
//   q0 = (1-delta_2/8)*q0 + (-q1*gx - q2*gy - q3*gz)*halfT;			// ������Ԫ����	 ��Ԫ��΢�ַ���	��Ԫ�������㷨�����ױϿ���
//   q1 = (1-delta_2/8)*q1 + (q0*gx + q2*gz - q3*gy)*halfT;
//   q2 = (1-delta_2/8)*q2 + (q0*gy - q1*gz + q3*gx)*halfT;
//   q3 = (1-delta_2/8)*q3 + (q0*gz + q1*gy - q2*gx)*halfT;			 
 												   
  q0 = q0 + (-q1*gx - q2*gy - q3*gz)*halfT;							// ������Ԫ����	 ��Ԫ��΢�ַ���	��Ԫ�������㷨��һ�����ⷨ
  q1 = q1 + (q0*gx + q2*gz - q3*gy)*halfT;
  q2 = q2 + (q0*gy - q1*gz + q3*gx)*halfT;         
  q3 = q3 + (q0*gz + q1*gy - q2*gx)*halfT;			   	   
//   mmx= q0 + (-q1*gx - q2*gy - q3*gz)*halfT;							// ������Ԫ����	 ��Ԫ��΢�ַ���	��Ԫ�������㷨��һ�����ⷨ
//   mmy= q1 + (q0*gx + q2*gz - q3*gy)*halfT;
//   mmz= q2 + (q0*gy - q1*gz + q3*gx)*halfT;         
//   mmf= q3 + (q0*gz + q1*gy - q2*gx)*halfT;		

  norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);		// ��������Ԫ
  q0 = q0 / norm;
  q1 = q1 / norm;
  q2 = q2 / norm;
  q3 = q3 / norm;
// printf(" %f; %f; %f; %f\n",q0,q1,q2,q3);
//   Q_ANGLE.Pitch = asin(-2*q1*q3 + 2*q0*q2)*57.3; 					//pitch
//   Q_ANGLE.Roll = atan2(2*q2*q3 + 2*q0*q1, -2*q1*q1 - 2*q2*q2+1)*57.3; 	// roll
//   Q_ANGLE.Yaw = -atan2(2*q1*q2 + 2*q0*q3, -2*q2*q2 - 2*q3*q3+1)*57.3; // yaw
//   //ת��Ϊŷ����,���ּ����������������м���������׼ȷ����
  Q_ANGLE.Roll =-asin(M13)*57.3; // pitch
	Q_ANGLE.Pitch = atan2(M23,M33)* 57.3; 	// roll
  Q_ANGLE.Yaw = -atan2(2*q1q2 + 2*q0q3,M11)* 57.3; // yaw
//   		printf("Q_ANGLE.Roll %f  Q_ANGLE.Pitch %f  AngEYaw % f\n",Q_ANGLE.Roll,Q_ANGLE.Pitch,Q_ANGLE.Yaw);

	
}




/****************************************************************************
* ��    �ƣ�Hcm_correction1()
* ��    �ܣ��شżƽ���
* ��ڲ�������
* ���ڲ�������
* ˵    ����
* ���÷������� 
****************************************************************************/ 


void  Hcm_correction()
{
	u8  hbuffer[8];
 	int16_t hcm_x, hcm_y, hcm_z;   				//��ȡ�شżƳ�ֵ
	
  IICreadBytes(HMC58X3_ADDR,HMC58X3_R_XM,6,hbuffer);
	hcm_x=(((int16_t)hbuffer[0]) << 8) | hbuffer[1];
	hcm_z=(((int16_t)hbuffer[2]) << 8) | hbuffer[3];
  hcm_y=(((int16_t)hbuffer[4]) << 8) | hbuffer[5];                      //�ֶ�����MATLAB��Բ�����㷨�������شżƵ���Χ�ų�ǿ��
	

	hcm_x= (hcm_x- HMC_OFFSET_X); 
 	hcm_y= (hcm_y- HMC_OFFSET_Y);
	hcm_z= (hcm_z- HMC_OFFSET_Z); //��ȥһ��ƫ��ֵ��Ȼ�����һ������������Y��Ϊ��׼ 
	
	Geomagnetic.X = ROTATE_XY[0]*hcm_x+ROTATE_XY[1]*hcm_y; 
 	Geomagnetic.Y = ROTATE_XY[2]*hcm_x+ROTATE_XY[3]*hcm_y;
	Geomagnetic.Z = ROTATE_YZ[2]*hcm_y+ROTATE_YZ[3]*hcm_z;
}
