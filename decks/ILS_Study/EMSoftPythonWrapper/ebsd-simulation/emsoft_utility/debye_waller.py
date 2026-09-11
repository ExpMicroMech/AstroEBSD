# debye Waller calculator
import pandas as pd
import numpy as np


def assign_debye(atoms: pd.DataFrame, dbw: float = None, temp: float = 293) -> pd.DataFrame:
    """This function assigns the debeye Waller factors based on either a input float for each atom, or based on an emprical model from this publication:
        The empirical debeye waller coefficients are read from a .txt file located in ./Support
    
    :param atoms: A pandas dataframe, containing the atomic species 
    :type atoms: pd.DataFrame
    :param dbw: Debye Waller Factor in the EM soft convention (nm)
    :type dbw: float
    :param temp: Temperature in K 
    :type temp: float"""
    df = df = pd.read_csv("./support/Debye_Waller.txt", sep=" ")
    df_el = df.drop_duplicates(subset="Element", keep="first")
    atoms_debye = pd.merge(
        atoms, df_el, left_on="element", right_on="Element", how="left"
    )
    # print(atoms_debye['Element'].isna())
    if not dbw == None:
        atoms_debye["dbw"] = dbw
    elif any(atoms_debye["Element"].isna()):
        print(
            """!!!!!!!!!!!!!!!!!! \n
              WARNING: Missing Element in Debye Waller Database \n
              !!!!!!!! \n"""
        )
    elif not "dbw" in atoms_debye:
        atoms_debye["dbw"] = calc_U(temp, atoms_debye)
    return atoms_debye


def calc_U(T : float, a_vals : float) -> list:
    """calc_U     
    Pass T in Kelvin and the parameters as described in the publication from Gao and Peng (1999)
    To fit to teh EM soft convention we need to divide by 100


    Args:
        T (_type_): _description_
        a_vals (_type_): _description_

    Returns:
        _type_: _description_
    """

    # df_work=df[df['Element']==element]
    # a=np.reshape(df_work[['a0','a1','a2','a3','a4']].values,5)
    dbw_list = []
    a_vals = a_vals.reset_index()
    for i in a_vals.index:
        a = np.reshape(a_vals.loc[i, ["a0", "a1", "a2", "a3", "a4"]].values, 5)
        dbw = (a[0] + a[1] * T + a[2] * (T**2) + a[3] * (T**3) + a[4] * (T**4)) / 100 # division by 100 to account for EM soft convention
        # dbw=0.001
        # dbw=0.0074
        dbw_list.append(dbw)
    # print(a)
    return dbw_list
