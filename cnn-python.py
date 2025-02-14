import xarray as xr
import torch
import captum
from torch.utils.data import DataLoader, random_split

import deep4downscaling.viz
import deep4downscaling.trans
import deep4downscaling.deep.loss
import deep4downscaling.deep.utils
import deep4downscaling.deep.models
import deep4downscaling.deep.train
import deep4downscaling.deep.pred
import deep4downscaling.metrics
import deep4downscaling.metrics_ccs
# Load predictor
predictor_filename = 'wa_pred.nc'
predictor = xr.open_dataset(predictor_filename)
# Load predictand
predictand_filename = '10km.nc'
predictand = xr.open_dataset(predictand_filename)
# Split data into training and test sets
years_train = ('1981', '2002')
years_test = ('2003', '2008')

x_train = predictor.sel(time=slice(*years_train))
y_train = predictand.sel(time=slice(*years_train))

x_test = predictor.sel(time=slice(*years_test))
y_test = predictand.sel(time=slice(*years_test))
x_train_stand = deep4downscaling.trans.standardize(data_ref=x_train, data=x_train)
y_mask = deep4downscaling.trans.compute_valid_mask(y_train) 
y_train_stack = y_train.stack(gridpoint=('lat', 'lon'))
y_mask_stack = y_mask.stack(gridpoint=('lat', 'lon'))

y_mask_stack_filt = y_mask_stack.where(y_mask_stack==1, drop=True)
y_train_stack_filt = y_train_stack.where(y_train_stack['gridpoint'] == y_mask_stack_filt['gridpoint'],
                                             drop=True)
                                             
y_train_nn = deep4downscaling.deep.utils.precipitation_NLL_trans(data=y_train_stack, threshold=0.9)
loss_function = deep4downscaling.deep.loss.NLLBerGammaLoss(ignore_nans=True)

x_train_stand_arr = deep4downscaling.trans.xarray_to_numpy(x_train_stand)
y_train_arr = deep4downscaling.trans.xarray_to_numpy(y_train_nn)

# Assuming your dataset is already created
train_dataset = deep4downscaling.deep.utils.StandardDataset(x=x_train_stand_arr, y=y_train_arr)

# Determine the total number of samples in the dataset
total_length = len(train_dataset)

# Calculate sizes for training and validation sets
train_size = int(0.9 * total_length)
valid_size = total_length - train_size

# Split the dataset into training and validation subsets
train_dataset, valid_dataset = random_split(train_dataset, [train_size, valid_size])

# Create DataLoaders
batch_size = 64
train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
valid_dataloader = DataLoader(valid_dataset, batch_size=batch_size, shuffle=True)

model_name = 'deepesd_pr'
model = deep4downscaling.deep.models.DeepESDpr(x_shape=x_train_stand_arr.shape,
                                               y_shape=y_train_arr.shape,
                                               filters_last_conv=1,
                                               stochastic=True)

num_epochs = 10000
patience_early_stopping = 20

learning_rate = 0.0001
optimizer = torch.optim.Adam(model.parameters(),
                             lr=learning_rate) 
                             
device = ('cuda' if torch.cuda.is_available() else 'cpu')

loss_function

DATA_PATH = './data/input'
FIGURES_PATH = './figures'
MODELS_PATH = './models'
ASYM_PATH = './data/asym'

train_loss, val_loss = deep4downscaling.deep.train.standard_training_loop(
                            model=model, model_name=model_name, model_path=MODELS_PATH,
                            device='cpu', num_epochs=num_epochs,
                            loss_function=loss_function, optimizer=optimizer,
                            train_data=train_dataloader, valid_data=valid_dataloader,
                            patience_early_stopping=patience_early_stopping)  
                            
                            
#Load the model weights into the DeepESD architecture
model.load_state_dict(torch.load(f'{MODELS_PATH}/{model_name}.pt'))

# Standardize
x_test_stand = deep4downscaling.trans.standardize(data_ref=x_train, data=x_test)

# Compute predictions
pred_test = deep4downscaling.deep.pred.compute_preds_ber_gamma(
                                x_data=x_test_stand, model=model,threshold=0.9,
                                device=device, var_target='precip',
                                mask=y_mask, batch_size=16) 
                                
# Visualize the predictions
deep4downscaling.viz.simple_map_plot(data=pred_test.mean('time'),
                                     colorbar='hot_r', var_to_plot='precip',
                                     output_path=f'./{FIGURES_PATH}/pred_mean.pdf') 
                                     
                                     
pred_test.to_netcdf('test.nc')                                                                                                                                      
