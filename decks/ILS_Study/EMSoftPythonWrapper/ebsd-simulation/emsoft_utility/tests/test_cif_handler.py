from emsoft_utility import cif_handler

def test_read_cif():
    lattice,atoms=cif_handler.read_cif('test_examples/mg.cif')

    assert lattice['gamma']==120
def test_crystal_system():
    outputs=[]
    outputs.append(cif_handler.assign_crystal_system(15))
    outputs.append(cif_handler.assign_crystal_system(41))
    outputs.append(cif_handler.assign_crystal_system(156))
    outputs.append(cif_handler.assign_crystal_system(195))
    assert outputs==[6,3,5,1]