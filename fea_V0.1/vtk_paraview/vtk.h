#pragma once
class vtk
{
private:
	bool dis_display;//��ʾ��ʵ�α��־λ TRUE����ʾ�α�  FALSE:����ʾ�α�
	float ratio;// ��ʾ�α�ʱ���α䱶��
public:
	vtk();
	vtk(bool dis_display, float ratio);
	void vtkoutput();//�ļ��������
};
