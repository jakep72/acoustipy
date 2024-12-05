from src.acoustipy.TMM import AcousticTMM
import torch
import matplotlib.pyplot as plt
import numpy as np
import time

def criterion(y1,y2,params):
    err = 1e12*torch.sum(torch.diff(y1-y2)**2)
    if params['vcl'] > params['tcl']:
        err = 2*err
    if params['fr'] > 1 or params['fr'] < 0:
        err = 2*err
    if params['phi'] > 1 or params['phi'] < 0.001:
        err = 2*err
    if params['tau'] > 1 or params['tau'] < 0.2:
        err = 2*err
    return(err)

class JCA(AcousticTMM):
    def __init__(self,
                 best_fr,
                 best_phi,
                 best_tau,
                 best_vcl,
                 best_tcl):

        super().__init__()
        self.fr = torch.nn.Parameter(best_fr)
        self.phi = torch.nn.Parameter(best_phi)
        self.tau = torch.nn.Parameter(best_tau)
        self.vcl = torch.nn.Parameter(best_vcl)
        self.tcl = torch.nn.Parameter(best_tcl)
        self.structure = AcousticTMM(incidence='Normal', air_temperature=20)
    
    def forward(self, thickness):
        layer = self.structure.Add_JCA_Layer(thickness, self.fr*1000000, self.phi, self.tau*5, self.vcl*500, self.tcl*500)
        tm = self.structure.assemble_structure(layer)
        a = self.structure.absorption(tm)[:,1].float()
        return a
    
    def string(self):
        return f'fr = {self.fr.item()*1000000} phi = {self.phi.item()} tau = {self.tau.item()*5} vcl = {self.vcl.item()*500} tcl = {self.tcl.item()*500}'
    
def quick_find(base_abs, thickness):
    fr = torch.linspace(10000,1000000,5)
    phi = torch.linspace(0.05,1,5)
    tau = torch.linspace(1,5,5)
    vcl = torch.linspace(10, 500, 5)
    tcl = torch.linspace(10, 500,5)
    best_err = 10
    for f in fr:
        for p in phi:
            for t in tau:
                for v in vcl:
                    for tc in tcl:
                        if tc >= v:
                            s = AcousticTMM(incidence='Normal',air_temperature=20)
                            l = s.Add_JCA_Layer(thickness,f,p,t,v,tc)
                            tm = s.assemble_structure(l)
                            a = s.absorption(tm)[:,1].float()
                            err = torch.sum(torch.diff(a-base_abs)**2)
                            
                            if err < best_err:
                                best_err = err
                                best_fr = f/1000000
                                best_phi = p
                                best_tau = t/5
                                best_vcl = v/500
                                best_tcl = tc/500

    return(best_fr, best_phi, best_tau, best_vcl, best_tcl)

def get_params(model):
    params = {}
    i=0
    for p in model.parameters():
        if i == 0:
            params['fr'] = p.item()
        if i == 1:
            params['phi'] = p.item()
        if i == 2:
            params['tau'] = p.item()
        if i == 3:
            params['vcl'] = p.item()
        if i == 4:
            params['tcl'] = p.item()
        i += 1
    return params

def clamper(model):
    i = 0
    for p in model.parameters():
        if i == 0:
            p.data.clamp_(0.00001, 1)
        if i==1:
            p.data.clamp_(0.001,0.999)
        elif i==2:
            p.data.clamp_(0.2,1)
        elif i==3:
            p.data.clamp_(0.002,1)
            # vcl = p.item()
        elif i==4:
            p.data.clamp_(.002,1)
            # tcl = p.item()
        i+=1

# Create an AcousticTMM object, specifying a diffuse sound field at 20C
structure = AcousticTMM(incidence='Normal',air_temperature=20)

# Define the layers of the material using various models
layer1 = structure.Add_JCA_Layer(25.4, 43000, .97, 1.2, 60, 134)


# Build the total transfer matrix of the structure + air gap
transfer_matrix = structure.assemble_structure(layer1)
a_const = structure.absorption(transfer_matrix)

# phi, tau, vcl, tcl = quick_find(a_const[:,1])
# Calculate the frequency dependent narrow band absorption coefficients
y = structure.absorption(transfer_matrix)[:,1].float()
# print(absorption)
s = time.time()
fr, phi, tau, vcl, tcl1 = quick_find(y, 25.4)
model = JCA(fr, phi, tau, vcl, tcl1)
# model = JCA(torch.randint(1,1000000,(1,))/1000000, torch.tensor(0.5), torch.randint(1,5,(1,))/5, torch.randint(1,500,(1,))/500, torch.randint(1,500,(1,))/500)

# y_pred = model(.0254)
# print(model.frequency)

# criterion = torch.nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
loss_mult = 1.0
loss_criteria = 1000

loss_curve = []
for t in range(100000):
    # Forward pass: Compute predicted y by passing x to the model
    y_pred = model.forward(25.4)
    
    # Compute and print loss
    loss = criterion(y_pred, y, get_params(model))
    if t == 0:
        init_loss = loss
    
    # print(loss)
    # if loss < prev_loss:
    #     best_params = model.string()
    #     print(t, loss.item(), best_params)
    
    # prev_loss = loss
    # if (loss / init_loss) > 20:
    #     init_loss = loss

    #     optimizer.param_groups[-1]['lr'] /= loss_mult
    #     loss_mult += 0.1
    #     loss_criteria = loss_criteria / 5
    #     print(g['lr'])
    #     print(loss_mult)
    # # print(init_loss)
    print(t, loss.item(), model.string())
    loss_curve.append(loss.item())
    if loss < 2000:
        for g in optimizer.param_groups:
            g['lr'] = 5e-4
    if loss < 250:
        for g in optimizer.param_groups:
            g['lr'] = 1e-4
    if loss < 10:
        for g in optimizer.param_groups:
            g['lr'] = 5e-5
    if loss < 8:
        break
    # if t % 100 == 0:
    
        # pred_layer = structure.Add_JCA_Layer(25.4, 105000,model.phi.item(), model.tau.item(), model.vcl.item()/(1e-6), model.tcl.item()/(1e-6))
        # tm = structure.assemble_structure(pred_layer)
        # a = structure.absorption(tm)
        # structure.plot_curve([a_const,a],['GT','Pred'])
        # plt.close()


    # Zero gradients, perform a backward pass, and update the weights.
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    # clamper(model)
    

    
e = time.time()
print(model.string())
print(e-s)








