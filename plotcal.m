[fn,fp,index] = uigetfile('*.csv');
if index == 0
    disp('No file selected')
else
    A = importdata(strcat(fp,fn));
    figure
    for i=2:size(A.data,2)
        f = A.data(:,i);
        bs = mean(f);
        df_scaled = (f-bs)/bs + 0.2*(i-2); % add 20% as a bias to plot all ROIs

        plot(A.data(:,1)/1000,df_scaled);
        hold on
    end
    legend(A.colheaders(1,2:end))
    legend boxoff
    title('Calcium signals')
    xlabel('time [s]')
    ylabel('\delta F / F')
end