## CORE BASE EXAMPLE 01

The will model a 600thk Slab of Grade 32MPa, 4.4m x 6.4m supported on a foundation material with a 
bearing capacity of 200kPa/mm. There will be tension only supports at the nodes:

| x mm        | y mm        |
| ----------- | ----------- |
| -1500.0     | -2500.0     |
|  1500.0     | -2500.0     |
| -1500.0     |  2500.0     |
|  1500.0     |  2500.0     |

![image](Core_Base_Example_01_mesh_plot.png)

### LOAD CASES

The load cases will be 'DL': Type 'DEAD', 'SDL': Type 'OTHER DEAD', 'LL': Type 'LIVE', 'LLRED': Type 'LIVE REDUCIBLE',
'EQX1': Type 'STATIC EQ X 1', 'EQY1': Type 'STATIC EQ Y 2'
