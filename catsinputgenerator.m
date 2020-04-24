function catsinputgenerator(hist)
%Generate the CATS probability and deviation files using a histogram of the dihedral coordinates as input. Each column of hist corresponds to a dihedral coordinate. Each line corresponds to a bin between -180  and 176.4

range=[-180:3.6:176.4]; %this range corresponds to the histogram bins.. the histogram program uses 100 bins at 3.6 degree increments starting from -180

options1=fitoptions('gauss1','Lower',[10 -200 0],'Upper',[max(hist(:)) 360 180]);
options2=fitoptions('gauss2','Lower',[10 -200 0 10 -200 0],'Upper',[max(hist(:)) 360 180 max(hist(:)) 360 180]);
options3=fitoptions('gauss3','Lower',[10 -200 0 10 -200 0 10 -200 0],'Upper',[max(hist(:)) 360 180 max(hist(:)) 360 180 max(hist(:)) 360 180]);

probfile = fopen('CATS_probabilities.out','w'); %open the file for writing
devfile = fopen('CATS_deviations.out','w');
redo='n';
i=1;
while i<=size(hist,2)
    disp(i)
    plot(range,hist(:,i))
    npeaks=input('npeaks\n')
	disp('choose start point')
    %textinput('Npeaks', npeaks)
    [x,y]=ginput(1);
    index=find(abs(range-x(1))<3.6);
    
    hist_sm=smoothdata(hist([index:length(range) 1:(index-1)],i),'gaussian',2); %smooth data using 1 neighboring point
    
    if(npeaks==1) %special cases for distribution/ranges
        
        f=fit([range(index:end) range(1:index-1)+360]',hist_sm,'gauss1',options1);
        plot(f,[range(index:end) range(1:index-1)+360],hist([index:length(range) 1:(index-1)],i))
        f
redo=input('redo?','s');
       if redo~='y'
            redo='n';
            fprintf(probfile,'%.4f, 1, 0, 0, 0, 0\n',f.b1);
            fprintf(devfile,'%.4f, 0, 0\n',f.c1);            
        end
        

    elseif npeaks==2
        
        
        f=fit([range(index:end) range(1:index-1)+360]',hist_sm,'gauss2',options2);
        
        p1=f.a1*f.c1/(f.a1*f.c1+f.a2*f.c2)
        p2=f.a2*f.c2/(f.a1*f.c1+f.a2*f.c2)
        plot(f,[range(index:end) range(1:index-1)+360],hist([index:length(range) 1:(index-1)],i))
        f
redo=input('redo?','s');
       if redo~='y'
            fprintf(probfile,'%.4f, %.4f, %.4f, %.4f, 0, 0\n',f.b1,p1,f.b2,p2);
        fprintf(devfile,'%.4f, %.4f, 0\n',f.c1,f.c2);
        end
        
        
    elseif npeaks==3
        
        
        f=fit([range(index:end) range(1:index-1)+360]',hist_sm,'gauss3',options3);
        
        f
        p1=f.a1*f.c1/(f.a1*f.c1+f.a2*f.c2 + f.a3*f.c3)
        p2=f.a2*f.c2/(f.a1*f.c1+f.a2*f.c2 + f.a3*f.c3)
        p3=f.a3*f.c3/(f.a1*f.c1+f.a2*f.c2 + f.a3*f.c3)
        plot(f,[range(index:end) range(1:index-1)+360],hist([index:length(range) 1:(index-1)],i))
       redo=input('redo?','s');
       if redo~='y'
         fprintf(probfile,'%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n',f.b1,p1,f.b2,p2,f.b3,p3);
        fprintf(devfile,'%.4f, %.4f, %.4f\n',f.c1,f.c2,f.c3);           
        end

        
    end
    
     if redo~='y'
       i=i+1;
    end
    
    close all
end
fclose(probfile);
fclose(devfile);
end

