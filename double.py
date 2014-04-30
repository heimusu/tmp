# -*- coding: utf-8 -*-
import cmath

# 定数
hbar = 1.054571628e-34
m0 = 9.10938215e-31
q = 1.602176487e-19
al = 0.063
m = q * al * m0
end_x = 6000.0
I = 0 + 1j

# 舞台設定など
delta_x = 10.0
x_0 = -5.0 * delta_x
E_0 = 0.155
S = 4500.0
L = 14.0

#無次元化パラメータ
Kr = 1 / (1.0e-9)
Er = ((hbar * hbar) * (Kr * Kr))/(2.0 * m)
k_0 = cmath.sqrt(E_0/Er)
V_0 = 0.3
d = 2.0
Sr = hbar/(q * Er)*1.0e15
s = S/Sr

#f(k)
def f(k):
    F = cmath.sqrt(cmath.sqrt(2.0 * delta_x * delta_x/cmath.pi)) * cmath.exp(-(k-k_0)*(k-k_0)*(delta_x)*(delta_x) - I*((k-k_0)*x_0))
    return F

def k2(k):
    res = cmath.sqrt(k*k - (V_0/Er))
    return res

def J(k):
    k_2 = k2(k)
    res = ((k_2 -  k) * cmath.exp(I * k * (L + d))) / (2.0 * k_2 * cmath.exp(-I * k_2 * (L + d)))
    return res

def H(k):
    k_2 = k2(k)
    res = ((k_2 + k) * cmath.exp(I * k * (L + d))) / (2.0 * k_2 * cmath.exp(I * k_2 * (L + d)));
    return res

def G(k):
    k_2 = k2(k)
    h = H(k)
    j = J(k)
    res = ((h * (k - k_2) * cmath.exp(I * k_2 * L)) / (2.0 * k * cmath.exp(-I * k * L))) + ((j * (k + k_2) * cmath.exp(-I * k_2 * L)) / (2.0 * k * cmath.exp(-I * k * L)))
    return res

def F(k):
    k_2 = k2(k)
    h = H(k)
    j = J(k)
    res = (h * (k + k_2)) * cmath.exp(I * k_2 * L) / (2.0 * k * cmath.exp(I * k * L)) + (j * (k - k_2)* cmath.exp(-I * k_2 * L)) / (2.0 * k * cmath.exp(I * k * L)) 
    return res



def D(k):
    k_2 = k2(k)
    f = F(k)
    g = G(k)
    res = ((f * (k_2 - k) * cmath.exp(I * k * d)) / (2.0 * k_2 * cmath.exp(-I * k_2 * d))) + ((g * (k_2 + k) * cmath.exp(-I * k * d)) / (2.0 * k_2 * cmath.exp(-I * k_2 * d)))
    return res



def C(k):
    k_2 = k2(k)
    f = F(k)
    g = G(k)
    res = ((f * (k_2 + k) * cmath.exp(I * k * d)) / (2.0 * k_2 * cmath.exp(I * k_2 * d))) + ((g * (k_2 - k) * cmath.exp(-I * k * d)) / (2.0 * k_2 * cmath.exp(I * k_2 * d)))
    return res


def B(k):
    k_2 = k2(k)
    c = C(k)
    d = D(k)
    res = (((k - k_2) * c) / (2.0 * k)) + (((k + k_2) * d) / (2.0 * k))
    return res


def A(k):
    k_2 = k2(k)
    c = C(k)
    d = D(k)
    res = (((k + k_2) * c) / (2.0 * k)) + (((k - k_2) * d) / (2.0 * k))
    return res


def R1(k):
    a = A(k)
    b = B(k)
    res = b/a
    return res


def A1(k):
    a = A(k)
    c = C(k)
    res = c/a
    return res

def B1(k):
    a = A(k)
    d = D(k)
    res = d/a
    return res


def T1(k):
    a = A(k)
    f = F(k)
    res = f/a
    return res

def R2(k):
    a = A(k)
    g = G(k)
    res = g/a
    return res

def A2(k):
    a = A(k)
    h = H(k)
    res = h/a
    return res

def B2(k):
    a = A(k)
    j = J(k)
    res = j/a
    return res

def T2(k):
    a = A(k)
    m = 1 + 0j
    res = m/a
    return res

if __name__ == "__main__":
    #std::ofstream ofs1("output.dat");
    #std::ofstream ofs2("fk.dat");
    #std::ofstream ofs3("tk.dat");
    #std::ofstream ofs4("etk.dat");
    #std::ofstream ofs5("dt.dat");
    f = open('output.dat','w')

    k1 = 0.25
    end_k = 0.8
    hx = 0.00625
    hk = 0.000125 #default:0.0005
    x = -1000.0
    res = 0 + 0j
    ex = int((end_x- x) / hx)
    ek = int((end_k - k1) / hk)
    k_sum = 0.0
    write1 = "k0 = "
    write2 = str(k_0)
    write3 = write1 + " " + write2

    #cout << "k0 = " << k_0 << endl
    f.write(write3)

    #set roop count
    N1 = int(-x / hx)
    N2 = int(d / hx)
    N3 = int((L- d) / hx)
    N4 = int(d / hx)
    N5 = int((end_x - (L + d)) / hx)

    tmp1 = abs(f(k1) * (cmath.exp(I * k1 * x) + R1(k1) * cmath.exp(-I * k1 * x)) * cmath.exp(-I * k1 * k1 * s))**2

    #入反射
    #for(nx = 0; nx <= N1; nx++){
    for nx in range(0,N1+1):
	k = k1 + hk
	res = 0 + 0j
	#for(jk = 1; jk < ek-1; jk+=2){
	for jk in range(1,ek-1,2):
	    res += 4.0 * f(k) * ( cmath.exp(I * k * x) + R1(k) * cmath.exp(-I * k * x) )* cmath.exp(-I * k * k * s)
	    k += hk * 2.0

	k = k1 + hk + hk;
	#for(jk = 2; jk < ek-1; jk+=2){
	for jk in range(2,ek-1,2):
	    res += 2.0 * f(k) * ( cmath.exp(I * k * x) + R1(k) * cmath.exp(-I * k * x) ) * cmath.exp(-I * k * k * s)
	    k += hk * 2.0

	res += f(k1) * ( cmath.exp(I * k1 * x) + R1(k1) * cmath.exp(-I * k1 * x) ) * cmath.exp(-I * k1 * k1* s)
	res += f(end_k) * ( cmath.exp(I * end_k * x) + R1(end_k) * cmath.exp(-I * end_k * x) ) * cmath.exp(-I * end_k * end_k * s)
	res = res * hk / 3.0
	res = res / sqrt(2.0 * cmath.pi)
	#ofs5 << setprecision(10) << x << " " << norm(res) << " " << real(res) << " " << imag(res) << endl;

	out = x + " " + str(abs(res)**2) + "\n"
	f.write(out)

	x += hx
	#台形則
	k_sum += ((tmp1 + (abs(res)**2) * hx) / 2.0)
	tmp1 = abs(res)**2

    #cout << "x " << x  << endl;
    #cout << "area " << k_sum << endl;
    print(x,str(x))
    print("area",str(k_sum))

    k_2 = k2(k1)
    tmp1 = abs((f(k1) * ( A1(k1) * cmath.exp(I * k_2 * x) +  B1(k1) * cmath.exp(-I * k_2 * x) ) * cmath.exp(-I * k1 * k1* s)))**2

    #障壁
    #for(nx = 0; nx <= N2; nx++){
    for nx in range(0,N2+1):
	k = k1 + hk
	res = 0+0j
	k_2 = k2(k)
	#for(jk = 1; jk<ek-1; jk+=2){
	for jk in range(1,ek-1,2):
	    res += 4.0 * f(k) * ( A1(k) * cmath.exp(I * k_2 * x) + B1(k) * cmath.exp(-I * k_2 * x) )* cmath.exp(-I * k * k * s)
	    k += hk * 2
	    k_2 = k2(k)

	k = k1 + hk + hk
	k_2 = k2(k)
	    #for(jk = 2; jk < ek-1; jk+=2){
	for jk in range(2,ek-1,2):
	    res += 2.0 * f(k) * ( A1(k) * cmath.exp(I * k_2 * x) + B1(k) * cmath.exp(-I * k_2 * x) ) * cmath.exp(-I * k * k * s)
	    k += hk * 2.0
	    k_2 = k2(k)
	
	k_2 = k2(k1)
	res += f(k1) * ( A1(k1) * cmath.exp(I * k_2 * x) + B1(k1) * cmath.exp(-I * k_2 * x) ) * cmath.exp(-I * k1 * k1* s)
	k_2 = k2(end_k)
	res += f(end_k) * ( A1(end_k) * cmath.exp(I * k_2 * x) + B1(end_k) * cmath.exp(-I * k_2 * x) ) * cmath.exp(-I * end_k * end_k * s)
	res  = res * hk / 3.0
	res = res / sqrt(2.0 * cmath.pi)
	#ofs5 << setprecision(10) << x << " " << norm(res) << " " << real(res) << " " << imag(res) << endl;
	out = x + " " + str(abs(res)**2) + "\n"
	f.write(out)
	x += hx

	k_sum += (tmp1 + norm(res)) * hx / 2.0
	tmp1 = abs(res)**2

	#cout << "x " << x << endl;
	#cout << "area " << k_sum << endl;
    print("x ",str(x))
    print("area ",str(k_sum))

    tmp1 = abs(f(k1) * ( T1(k1) * cmath.exp(I * k1 * x) + R2(k1) * cmath.exp(-I * k1 * x) ) * cmath.exp( - I * k1 * k1* s) )**2
    #//入反射2
    #for(nx = 0; nx <= N3; nx++){
    for nx in range(0,N3+1):
	k = k1 + hk
	res = 0 + 0j
	#for(jk = 1; jk<ek-1; jk+=2){
	for jk in range(1,ek-1,2):
	    res += 4.0 * f(k) * ( T1(k) * cmath.exp(I * k * x) + R2(k) * cmath.exp( -I * k * x ) ) * cmath.exp( - I * k * k * s)
	    k += hk * 2.0
	
	k = k1 + hk + hk
	#for(jk = 2; jk < ek-1; jk+=2){
	for jk in range(2,ek-1,2):
	    res += 2.0 * f(k) * ( T1(k) * cmath.exp(I * k * x) + R2(k) * cmath.exp(-I * k * x) ) * cmath.exp( - I * k * k * s)
	    k += hk * 2.0

	res += f(k1) * ( T1(k1) * cmath.exp(I * k1 * x) + R2(k1) * cmath.exp(-I * k1 * x )) * cmath.exp( - I * k1 * k1* s)
	res += f(end_k) * ( T1(end_k)  * cmath.exp(I * end_k * x) + R2(end_k) * cmath.exp(-I * end_k * x)) * cmath.exp( - I * end_k * end_k * s)
	res  = res * hk / 3.0
	res = res / sqrt(2.0 * cmath.pi)
	#ofs5 << setprecision(10) << x << " " << norm(res) << " " << real(res) << " " << imag(res) << endl;
	out = x + " " + str(abs(res)**2) + "\n"
	f.write(out)
	x += hx

	k_sum += (tmp1 + abs(res)**2) * hx / 2.0
	tmp1 = abs(res)**2

    #cout << "x " << x << endl;
    #cout << "area " << k_sum << endl;
    print("x ",str(x))
    print("area ",str(k_sum))

    k_2 = k2(k1)
    tmp1 = abs(f(k1) * (A2(k1) * cmath.exp(I * k_2 * x) + B2(k1) * cmath.exp(-I * k_2 * x)) * cmath.exp(-I * k1 * k1* s))**2

    #//障壁2
    #for(nx = 0; nx <= N4; nx++){
    for nx in range(0,N4+1):
	k = k1 + hk
	res = 0 + 0j 
	k_2 = k2(k)
	#for(jk = 1; jk<ek-1; jk+=2){
	for jk in range(1,ek-1,2):
	    res += 4.0 * f(k) * (A2(k) * cmath.exp(I * k_2 * x) + B2(k) * cmath.exp(-I * k_2 * x))* cmath.exp(-I * k * k * s)
	    k += hk * 2.0
	    k_2 = k2(k)

	k = k1 + hk + hk
	k_2 = k2(k)
	#for(jk = 2; jk < ek-1; jk+=2){
	for jk in range(2,ek-1,2):
	    res += 2.0 * f(k) * (A2(k) * cmath.exp(I * k_2 * x) + B2(k) * cmath.exp(-I * k_2 * x)) * cmath.exp(-I * k * k * s)
	    k += hk * 2.0
	    k_2 = k2(k)

	k_2 = k2(k1)
	res += f(k1) * (A2(k1) * cmath.exp(I * k_2 * x) + B2(k1) * cmath.exp(-I * k_2 * x)) * cmath.exp(-I * k1 * k1* s)
	k_2 = k2(end_k)
	res += f(end_k) * (A2(end_k) * cmath.exp(I * k_2 * x) + B2(end_k) * cmath.exp(-I * k_2 * x)) * cmath.exp(-I * end_k * end_k * s)
	res  = res * hk / 3.0
	res = res / sqrt(2.0 * cmath.pi)
	#ofs5 << setprecision(10) << x << " " << norm(res) << " " << real(res) << " " << imag(res) << endl;
	out = x + " " + str(abs(res)**2) + "\n"
	f.write(out)
	x += hx
	k_sum += (tmp1 +abs(res)**2) * hx / 2.0
	tmp1 = abs(res)**2


    #cout << "x " << x  << endl;
    #cout << "area " << k_sum << endl;
    print("x ",str(x))
    print("area ",str(k_sum))

    tmp1 = abs(f(k1) * T2(k1) * cmath.exp(I * k1 * x - I * k1 * k1* s))**2
    #//透過
    #for(nx = 0; nx <= N5; nx++){
    for nx in range(0,N5+1):
	k = k1 + hk
	res = 0 + 0j 
	#for(jk = 1; jk<ek-1; jk+=2){
	for jk in range(1,ek-1,2):
	    res += 4.0 * f(k) * T2(k) * cmath.exp(I * k * x - I * k * k * s)
	    k += hk * 2.0

	k = k1 + hk + hk
	#for(jk = 2; jk < ek-1; jk+=2){
	for jk in range(2,ek-1,2):
	    res += 2.0 * f(k) * T2(k) * cmath.exp(I * k * x - I * k * k * s)
	    k += hk * 2.0
	
	res += f(k1) * T2(k1) * cmath.exp(I * k1 * x - I * k1 * k1* s)
	res += f(end_k) * T2(end_k)  * cmath.exp(I * end_k * x - I * end_k * end_k * s)
	res  = res * hk / 3.0
	res = res / sqrt(2.0 * cmath.pi)
	out = x + " " + str(abs(res)**2) + "\n"
	f.write(out)
	x += hx
	k_sum += (tmp1 +abs(res)**2) * hx / 2.0
	tmp1 = abs(res)**2
    
    #cout << "x " << x  << endl;
    #cout << "area " << k_sum << endl;
    print("x ",str(x))
    print("area ",str(k_sum))

    f.close()
