#ifndef mathfem_H_
#define mathfem_H_
#include "DoubleMatrix.h"
//���أ����ֵС��0������-1�����򷵻�1
inline double sgn(double i)
{
	return (i < 0. ? -1. : 1.);
}

//���ֵ
inline int max_fem(const int a, const int b) { return a >= b ? a : b; }
inline double max_fem(const double a, const double b) { return a >= b ? a : b; }

//��Сֵ
inline int min_fem(const int a, const int b) { return a <= b ? a : b; }
inline double min_fem(const double a, const double b) { return a <= b ? a : b; }

/*
 *������η���(a*x^3+b*x^2+c*x+d=0)
 *@r1:��һ��
 *@r2:�ڶ���
 *@r3:������
 *@nroot:������Ŀ
 */
void cubic(double a, double b, double c, double d, double *r1, double *r2, double *r3, int *nroot);
/**
 * \brief ���ԳƷ���2*2 or 3*3��������ֵ��������������Jacobi�������㾫ȷ���ٶȽ���
 * \param vec ��������
 * \param val ����ֵ
 * \param src ����ĶԳƾ���(src���ڿ������Ƿ�����������)
 */
void eigOfSymmtricMatrix(DoubleMatrix& vec, DoubleArray& val, DoubleMatrix src);

/**
 * \brief �ж�ƽ�����߶�1���߶�2�Ƿ��ཻ
 * \param x1 �߶�1��һ��ĺ�����
 * \param y1 �߶�1��һ���������
 * \param x2 �߶�1�ڶ���ĺ�����
 * \param y2 �߶�1�ڶ����������
 * \param x3 �߶�2��һ��ĺ�����
 * \param y3 �߶�2��һ���������
 * \param x4 �߶�2�ڶ���ĺ�����
 * \param y4 �߶�2�ڶ����������
 * \return �ཻ����true,���뽻����false
 */
bool judgeInsection(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

/*Definitions of useful mathematical constants
* M_E - e
* M_LOG2E - log2(e)
* M_LOG10E - log10(e)
* M_LN2 - ln(2)
* M_LN10 - ln(10)
* M_PI - pi
* M_PI_2 - pi / 2
* M_PI_4 - pi / 4
* M_1_PI - 1 / pi
* M_2_PI - 2 / pi
* M_2_SQRTPI - 2 / sqrt(pi)
* M_SQRT2 - sqrt(2)
* M_SQRT1_2 - 1 / sqrt(2)
*/

#define M_E        2.71828182845904523536
#define M_LOG2E    1.44269504088896340736
#define M_LOG10E   0.434294481903251827651
#define M_LN2      0.693147180559945309417
#define M_LN10     2.30258509299404568402
#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#define M_PI_4     0.785398163397448309616
#define M_1_PI     0.318309886183790671538
#define M_2_PI     0.636619772367581343076
#define M_2_SQRTPI 1.12837916709551257390
#define M_SQRT2    1.41421356237309504880
#define M_SQRT1_2  0.707106781186547524401
#endif