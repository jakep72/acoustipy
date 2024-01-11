# Saving layers to a database
Save the parameters of individual layers to a built-in SQLite database.

```python
from acoustipy import AcousticTMM, AcoustiBase


structure = AcousticTMM(incidence='Normal',air_temperature = 20)

# For each Add_XXX_Layer method, enable the save_layer parameter and give it a unique layer name
db = structure.Add_DB_Layer(25.4, 50000,save_layer=True,layer_name='test_DB')
dbm = structure.Add_DBM_Layer(25.4, 50000,save_layer=True,layer_name='test_DBM')
jca = structure.Add_JCA_Layer(25.4, 50000,.90,1.2,80,110,save_layer=True,layer_name='test_JCA')
jcal = structure.Add_JCAL_Layer(25.4, 50000,.90,1.2,80,110,50,save_layer=True,layer_name='test_JCAL')
jcapl = structure.Add_JCAPL_Layer(25.4, 50000, .90, 1.2, 80, 110, 50, 70, 65,save_layer=True,layer_name='test_JCAPL')
horoshenkov = structure.Add_Horoshenkov_Layer(40, .90, 147, .325, save_layer=True, layer_name='test_horoshenkov')
biot_limp = structure.Add_Biot_Limp_Layer('JCA',25.4, 50000,200, .90, 1.2, 80, 110,save_layer=True,layer_name='test_biot_limp')
biot_rigid = structure.Add_Biot_Rigid_Layer('JCA',25.4, 50000,200, .90, 1.2, 80, 110,save_layer=True,layer_name='test2_biot_rigid')
screen = structure.Add_Resistive_Screen(2.54, 50000, .90,save_layer=True,layer_name='test_screen')
maa = structure.Add_MAA_MPP_Layer(2.54, .5, 1,save_layer=True,layer_name='test_maa_mpp')
ef_mpp = structure.Add_MPP_EF_Layer(2.54, .5, 1,save_layer=True,layer_name='test_ef_mpp')

# Pull the layer information from the database and save to a .csv file
s = AcoustiBase()
s.summarize_layers()
```