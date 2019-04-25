clear all;

%D regulator
%classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, ManVal) 
%D regulator
p=classPID(0, 0, 1, 30, 1, 100, -100, 1, 1)

%SetAutoMan(obj, AutoMan, ManVal)
p.SetAutoMan(1,1)

%PID
%reTune(obj, K, Ti, Kd, Td)
p.reTune(1, 30, 0.5, 30)

stpt = 1;
pv=0;
u=0;

%obiekt = tf(1,[30 1]);
for i=1:1:10
  u =  p.calc(pv,stpt)
  %wyliczenie wyjscia obiektu na nastepny krok wykorzystujac u
  
end
