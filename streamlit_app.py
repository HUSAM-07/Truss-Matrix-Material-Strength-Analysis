import streamlit as st
import math
import numpy as np

def truss_analysis(tn, te, xco, yco, A, E, elements, supports, loads):
    np.set_printoptions(3, suppress=True)

    snofel = []
    enofel = []
    lenofel = []
    elcon = []
    cosofel = []
    sinofel = []

    for i in range(te):
        a, b = elements[i]
        x1, y1 = xco[a-1], yco[a-1]
        x2, y2 = xco[b-1], yco[b-1]
        l = math.sqrt((x2-x1)**2 + (y2-y1)**2)

        if l == 0:
            st.warning(f"Element {i+1} has zero length. Please adjust node positions.")
            return None, None, None, None, None, None, None, None

        con = A*E/l
        cos = (x2-x1)/l
        sin = (y2-y1)/l

        snofel.append(a)
        enofel.append(b)
        lenofel.append(l)
        elcon.append(con)
        cosofel.append(cos)
        sinofel.append(sin)

    elstmat = []

    for i in range(te):
        cc = float(cosofel[i])**2
        ss = float(sinofel[i])**2
        cs = float(cosofel[i])*float(sinofel[i])

        mat = elcon[i]*np.array([[cc, cs, -cc, -cs],
                                [cs, ss, -cs, -ss],
                                [-cc, -cs, cc, cs],
                                [-cs, -ss, cs, ss]])

        elstmat.append(mat)

    gstmatmap = []

    for i in range(te):
        m = snofel[i]*2
        n = enofel[i]*2
        add = [m-1, m, n-1, n]
        gmat = np.zeros((tn*2, tn*2))
        elmat = elstmat[i]
        for j in range(4):
            for k in range(4):
                a = add[j]-1
                b = add[k]-1
                gmat[a, b] = elmat[j, k]
        gstmatmap.append(gmat)

    GSM = np.zeros((tn*2, tn*2))
    for mat in gstmatmap:
        GSM = GSM + mat

    displist = []
    forcelist = []
    for i in range(tn):
        a = str('u')+str(i+1)
        displist.append(a)
        b = str('v')+str(i+1)
        displist.append(b)
        c = str('fx')+str(i+1)
        forcelist.append(c)
        d = str('fy')+str(i+1)
        forcelist.append(d)

    dispmat = np.ones((tn*2, 1))
    for support in supports:
        supn, condition = support
        if condition in ['P', 'p']:
            dispmat[supn*2-2, 0] = 0
            dispmat[supn*2-1, 0] = 0
        elif condition in ['H', 'h']:
            dispmat[supn*2-2, 0] = 0
        elif condition in ['V', 'v']:
            dispmat[supn*2-1, 0] = 0
        else:
            st.warning('Invalid support condition')

    for load in loads:
        lon, fx, fy = load
        forcelist[lon*2-2] = fx
        forcelist[lon*2-1] = fy

    dispresult = np.matmul(np.linalg.inv(GSM), forcelist)
    rin = 0
    for i in range(tn*2):
        if dispmat[i, 0] == 1:
            dispmat[i, 0] = dispresult[rin, 0]
            rin = rin + 1

    rcdlist = []
    for i in range(tn*2):
        if dispmat[i, 0] == 0:
            rcdlist.append(i)

    rrgsm = np.delete(GSM, rcdlist, 0)
    crgsm = np.delete(rrgsm, rcdlist, 1)
    rgsm = crgsm
    rforcemat = np.delete(forcelist, rcdlist, 0)

    dispresult = np.matmul(np.linalg.inv(rgsm), rforcemat)

    forcelist = np.matmul(GSM, dispmat)

    newxco = []
    newyco = []
    count = 0
    for i in range(tn):
        k = xco[i]+dispmat[count, 0]
        newxco.append(k)
        count = count + 1
        l = yco[i]+dispmat[count, 0]
        newyco.append(l)
        count = count + 1

    newlenofel = []
    for i in range(te):
        a, b = snofel[i], enofel[i]
        x1 = float(newxco[a-1])
        y1 = float(newyco[a-1])
        x2 = float(newxco[b-1])
        y2 = float(newyco[b-1])
        l = math.sqrt((x2-x1)**2 + (y2-y1)**2)
        newlenofel.append(l)

    elstrain = np.zeros((te, 1))
    for i in range(te):
        elstrain[i, 0] = (newlenofel[i]-lenofel[i])/(lenofel[i])

    elstress = np.zeros((te, 1))
    for i in range(te):
        elstress[i, 0] = E * elstrain[i, 0]

    eforce = np.zeros((te, 1))
    for i in range(te):
        eforce[i, 0] = A * elstress[i, 0]

    return GSM, dispmat, forcelist, newxco, newyco, elstrain, elstress, eforce

def main():
    # Navigation
    pages = ["Truss Matrix Info", "Run Analysis"]
    page = st.sidebar.selectbox("Select a page:", pages)

    

    if page == "Truss Matrix Info":
        st.write("# Welcome to Truss Analysis App")
        st.write("This app helps you analyze truss structures and provides information about the truss matrix.")
        st.write("# Truss Matrix Information")
        st.write("A truss is a structure that consists of members organized to form triangles. The analysis of truss structures involves solving for displacements and forces using the stiffness matrix.")
        

        st.subheader("Global Stiffness Matrix (GSM):")
        st.write("The GSM is a square matrix that represents the stiffness characteristics of the entire truss structure. It is assembled by combining the stiffness matrices of individual truss elements.")

        st.subheader("Displacement Matrix:")
        st.write("The displacement matrix contains the displacements of each degree of freedom at every node. It is calculated by solving the equilibrium equations based on boundary conditions and loads.")

        st.subheader("Force Matrix:")
        st.write("The force matrix represents the forces applied to each degree of freedom at every node. It is determined by applying external loads to the truss structure.")

        st.warning("This app is made by Krithika Rathi", icon='😊')

    elif page == "Run Analysis":

        tn = st.number_input('Enter the total number of nodes:', value=3, min_value=1, key='tn')
        te = st.number_input('Enter the total number of Elements:', value=3, min_value=1, key='te')

        xco = []
        yco = []
        for i in range(tn):
            x = st.number_input(f'Enter the x co-ordinate of node {i+1} in mm:', key=f'x_{i}', value=0.0)
            y = st.number_input(f'Enter the y co-ordinate of node {i+1} in mm:', key=f'y_{i}', value=0.0)
            xco.append(x)
            yco.append(y)

        A = st.number_input('Enter the Area of cross section in mm2:', value=10.0, key='A')
        E = st.number_input('Enter the Modulus of Elasticity in N/mm2:', value=100.0, key='E')

        elements = []
        for i in range(te):
            a = st.number_input(f'Enter the Start node of element {i+1}:', key=f'start_{i}', min_value=1, max_value=tn)
            b = st.number_input(f'Enter the End node of element {i+1}:', key=f'end_{i}', min_value=1, max_value=tn)
            elements.append((a, b))

        supports = []
        loads = []

        tsupn = st.number_input('Enter the total number of nodes having supports:', min_value=0, max_value=tn, key='tsupn')
        supcondition = ['P = pinned',
                        'H = Horizonal restrained (vertical is free to move)',
                        'V = Vertical restrained (Horizontal is free to move)']

        for i in range(tsupn):
            supn = st.number_input('Enter the node number of support:', min_value=1, max_value=tn, key=f'supn_{i}')
            for a in supcondition:
                st.write(a)
            condition = st.text_input('Enter the condition of the support:', key=f'condition_{i}')
            supports.append((supn, condition))

        tlon = st.number_input('Enter the total number of loaded nodes:', min_value=0, max_value=tn, key='tlon')
        for i in range(tlon):
            lon = st.number_input('Enter the node number of Loading:', min_value=1, max_value=tn, key=f'lon_{i}')
            fx = st.number_input('Enter the Horizontal load at this node in N:', key=f'fx_{lon}', value=0.0)
            fy = st.number_input('Enter the Vertical load at this node in N:', key=f'fy_{lon}', value=0.0)
            loads.append((lon, fx, fy))

        if st.button("Run Analysis"):
            result = truss_analysis(tn, te, xco, yco, A, E, elements, supports, loads)

            if result is not None:
                GSM, dispmat, forceresult, newxco, newyco, elstrain, elstress, eforce = result

                st.subheader('Results:')
                st.write('\n\nGlobal Stiffness Matrix of the Truss')
                st.write(GSM)

                st.write('\n\nDisplacement matrix of nodes')
                st.write(dispmat)

                st.write('\n\nForce matrix of nodes')
                st.write(forceresult)

                st.write('\n\nNew coordinates of nodes')
                st.write(list(zip(newxco, newyco)))

                st.write('\n\nStrain in the elements')
                st.write(elstrain)

                st.write('\n\nStress in the elements')
                st.write(elstress)

                st.write('\n\nForce in the elements')
                st.write(eforce)

if __name__ == "__main__":
    main()
