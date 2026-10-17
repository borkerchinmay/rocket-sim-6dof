function C = quat2DCM(q)
%QUAT2DCM  Compute Direction Cosine Matrix from unit quaternion.
%
%  Input:  q = [q0; q1; q2; q3]  scalar-first convention
%          q represents rotation from BODY to NED
%          C = C_BN  transforms vectors from NED → Body
%
%  C_BN * v_NED = v_body
%  C_NB = C_BN'   transforms body → NED

q0 = q(1);  q1 = q(2);  q2 = q(3);  q3 = q(4);

% Pre-compute products for efficiency
q0q0 = q0*q0;  q1q1 = q1*q1;  q2q2 = q2*q2;  q3q3 = q3*q3;
q0q1 = q0*q1;  q0q2 = q0*q2;  q0q3 = q0*q3;
q1q2 = q1*q2;  q1q3 = q1*q3;  q2q3 = q2*q3;

C = [q0q0+q1q1-q2q2-q3q3,   2*(q1q2+q0q3),         2*(q1q3-q0q2);
     2*(q1q2-q0q3),           q0q0-q1q1+q2q2-q3q3,   2*(q2q3+q0q1);
     2*(q1q3+q0q2),           2*(q2q3-q0q1),         q0q0-q1q1-q2q2+q3q3];
end
