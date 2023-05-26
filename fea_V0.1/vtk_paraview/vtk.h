#pragma once
class vtk
{
private:
	bool dis_display;//显示真实形变标志位 TRUE：显示形变  FALSE:不显示形变
	float ratio;// 显示形变时，形变倍数
public:
	vtk();
	vtk(bool dis_display, float ratio);
	void vtkoutput();//文件输出函数
};
