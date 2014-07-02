%% Main File for Running Halo Assay Program
clear all

if exist('gui', 'var')
    try
        delete(gui)
    catch
    end
    clear gui
end
gui = Initialization;