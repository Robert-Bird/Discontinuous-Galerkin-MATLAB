function ImposedDisplacements = ImposedDisplacements(X,Y)
%IMPOSEDDISPLACEMENTS
%    IMPOSEDDISPLACEMENTS = IMPOSEDDISPLACEMENTS(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    28-Mar-2022 16:43:01

t2 = X.*pi;
t3 = Y.*pi;
t4 = sin(t2);
t5 = sin(t3);
t6 = t4.*t5;
ImposedDisplacements = [t6;t6];
