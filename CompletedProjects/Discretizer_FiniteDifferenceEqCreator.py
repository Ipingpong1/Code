# Script which returns the finite difference approximation of any given equation
# 2021

eq = 'dh/dt = - du/dx - dv/dy'
eqarr = eq.split(' ')
func = ['h[t, x, y]', 'u[t, x, y]', 'v[t, x, y]']

def diffToFDM(str, fun, wr, r):
    funp = fun[0:fun.index(wr)+1] + "+1" + fun[fun.index(wr)+1:-1]+')'
    funm = fun[0:fun.index(wr)+1] + "-1" + fun[fun.index(wr)+1:-1]+')'

    if(r == 1):
        return ('('+funp + ' - ' + fun + ')/d' + wr)
    if(r == 2):
        return ('('+funp+' + '+'2 * ' + fun + ' - ' + funm +')/d'+wr+"^2")

def diffToBDM(str, fun, wr, r):
    funp = fun[0:fun.index(wr)+1] + "+1" + fun[fun.index(wr)+1:-1]+')'
    funm = fun[0:fun.index(wr)+1] + "-1" + fun[fun.index(wr)+1:-1]+')'

    if(r == 1):
        return ('('+fun + ' - ' + funm + ')/d' + wr)
    if(r == 2):
        return ('('+funp+' + '+'2 * ' + fun + ' - ' + funm +')/d'+wr+"^2")

def isdiff(str):
    for fun in func:
        if(str.find('d'+fun[0]+'/d')!=-1):
            return diffToFDM(str, fun, str[-1], 1)

        if(str.find('d^2'+fun[0]+'/d')!=-1):
            return diffToFDM(str, fun, str[-3], 2)

    return str

final = ''

for i in range(0, len(eqarr)):
    eqarr[i] = isdiff(eqarr[i])
    final+=eqarr[i]+' '

print(final)


