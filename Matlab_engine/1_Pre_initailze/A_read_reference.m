function [Ureff,Treff,Roreff,mureff,c] = A_read_reference
    A = importdata('input/Ref.dat',' ',1);
    a = A.data;

    Ureff  = a(1);         %�ο��ٶ�
    Treff  = a(2);         %�ο��¶�
    Roreff = a(3);         %�ο��ܶ�
    mureff = a(4);         %�ο�ճ��ϵ��
    c      = a(5);        %�ο��ҳ�
end