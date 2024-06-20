import modin.pandas as pandas
import ray


ray.init(_plasma_directory="/tmp") # setting to disable out of core in Ray

# pip install python-calamine
horse_residual_height_data=pandas.read_excel("test_env/MapResults_rHeight.xls", index_col=None, engine="calamine")
print(horse_residual_height_data.head())

horse_residual_height_data['mlog10_pSNP8'] = horse_residual_height_data['mlog10_pSNP8'].fillna(9999)

horse_residual_height_data = horse_residual_height_data.sort_values(by=['mlog10_pSNP8'],ascending=True)
print("Head")
print(horse_residual_height_data.head(100))

pos_top=50
pos_bot=70
print(f'Viewing positions {pos_top} - {pos_bot}')
print(horse_residual_height_data.iloc[pos_top:pos_bot,])

print("Tail")
horse_residual_height_data = horse_residual_height_data.sort_values(by=['mlog10_pSNP8'],ascending=False)
print(horse_residual_height_data.head(100))

