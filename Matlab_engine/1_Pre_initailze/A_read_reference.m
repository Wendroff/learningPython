function [Ureff,Treff,Roreff,mureff,c] = A_read_reference
    A = importdata('input/Ref.dat',' ',1);
    a = A.data;

    Ureff  = a(1);         %参考速度
    Treff  = a(2);         %参考温度
    Roreff = a(3);         %参考密度
    mureff = a(4);         %参考粘性系数
    c      = a(5);        %参考弦长
end