function po = trilateration3D(p,r)
% p -> 3D position of N anchors
% r -> distance measurement from tag to the N anchors
N = length(r);
sum_a = zeros(3,1);
sum_B = zeros(3,3);
sum_c = zeros(3,1);
sum_H = zeros(3,3);
sum_qtq1 = 0;
sum_qtq2 = 0;
for i = 1:N
    pi = p(:,i);
    ri = r(i);
    sum_a = sum_a + pi*pi'*pi - ri^2*pi;
    sum_B = sum_B -2*pi*pi' - pi'*pi*eye(3) + ri^2*eye(3);
    sum_c = sum_c + pi;
    sum_H = sum_H + pi*pi';
    sum_qtq1 = sum_qtq1 + pi'*pi;
    sum_qtq2 = sum_qtq2 + ri^2;
end

a = (1/N)*sum_a;
B = (1/N)*sum_B;
c = (1/N)*sum_c;
f = a + B*c + 2*c*c'*c;
H = (-2/N)*sum_H + 2*c*c';
q = -H^(-1)*f;
po = q + c;
end