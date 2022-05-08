# Packages:
from numpy import inf, nan, sin, cos, pi, sqrt, diff, std, diag, argmin
from numpy import mean, median, nanargmax, zeros, ones, ceil, delete, shape
from numpy import arange, array, size, vstack, copy, loadtxt, where, savetxt, linspace
from numpy import min, max, sum, float, round, int
import math, sys, time
import numpy as np
import pyfits, pylab, scipy 
import matplotlib.pyplot as plt
import scipy.ndimage as snd
from matplotlib.colors import LogNorm
from scipy.misc import imsave
# Functions:


def Lucky_Imaging(data, FSR, plot=None, save=None): 
    """------------------------------------------- FUNCTION ------------------------------------------------:
    This function takes a cube of very short (close to the defraction limit) exposure images an sort them by 
    quality and combine them to a single image.  
    ----------INPUT:
    data           : Cube of short exposures.
    FSR            : Frame Selection Rate (FSR); procentages of images that is used for the final image.
    ---------OUTPUT:
    Q              : The quality """
    print '------------------------------------------------------------------------- Lucky_Imaging'
    start_time = time.time()  # Take time

    r     = 10
    r_min = 10
    r_max = 13
    r_opt = 5

    x_star = 300; y_star = 300
    
    print ('Filter done in time: %s s' % (time.time() - start_time))
    return Q

# %% Optimal radii:
# % Find bedste radius, kald denne vaerdi r_op
# % radius = optimal_radius(img,px,py,rInn,rOut,pL_);


# % AperturePhotometry(data(:,:,8),398,385,4,10)


# center      = Find_Higest_Pixel(image,y_star,x_star,d);
# % center3     = Find_Only_Higest_Pixel(data,yStar,xStar,d);
# data_shift  = Shift_Frames(image,center);
# % data_shift3 = Shift_Frames(data,center3);
# i1          = Quality_Aperture(image,y_star,x_star,center,d);
# % i2          = Quality_Gauss_Fit(data,yStar,xStar,center,d);
# % i3          = Quality_Aperture(data,yStar,xStar,center3,d);
# Plots(image,data_shift,i1);

# toc 

# end

# function radius = optimal_radius(image,x_star,y_star,r_in,r_out)

# r_1      = 1;
# dr       = 0.5;
# r_2      = 12;
# r_opti   = [];
# index    = [];

#   for i = 1:length(x_star)      % Find optimal radius for all stars
      
#         SNR_rad     = [];
#         f_radius    = [];
        
#     for r=r_1:dr:r_2
        
#          [n_pix N_sky f] = aperture(image,x_star(i),y_star(i),r,r_in,r_out,pL_);
#          SNR_rad = [SNR_rad signal_to_noise(f,n_pix,N_sky)];
#          f_radius = [f_radius f];

#     end
    
#     if sum(gradient(SNR_rad) < 0) > 5 && sum(gradient(SNR_rad) > 0) > 5
        
#         [snrM I] = max(SNR_rad); 
#         index = [index; i];
#         r_opti = [r_opti; I * dr];
# %       figure % plot flux vs radius
# %       plot(r_1:dr:r_2,f_radius,'.')
# %       xlabel('radius'); ylabel('flux [ ADU ]');
# %       figure
# %       plot(r_1:dr:r_2,SNR_rad,'.') % plus SNR vs radius
# %       xlabel('radius'); ylabel('SNR');
 
#     end
  
#   end
  
#   radius = 3*median(r_opti)-2*mean(r_opti);
  
# end
# function [center] = Find_Only_Higest_Pixel(data,yStar,xStar,d)
# % This function takes the a data cube and finds the center for a reference star
# % in each frame. The function takes the raw cube (each frame is a fits
# % file), the position of the reference star and a distance to the square.

# center_x = [];
# center_y = [];
# center   = [];
# for i = 1:length(data(1,1,:));
    
# square          = data(yStar-d:yStar+d,xStar-d:xStar+d,i); % Udsnit i hvert billede

# [~,x_max_i]     = max(max(square));         % Søjlen og værdien for max pixel
# [~,y_max_i]     = max(square(:,x_max_i));   % Rækken for max pixel                   

# center_y_i      = y_max_i;                  % Søjlerne samles i vektor  
# center_x_i      = x_max_i;                  % Rækkerne samles i vektor

# center_x        = [center_x center_x_i];
# center_y        = [center_y center_y_i];
# center          = [center_y' center_x'];    % center opstilles i 2 søjler

# % figure(2)
# % imshow(square,[0,10000])
# % hold on
# % plot(center_x_i,center_y_i,'b*')
# % set(gca,'box','on',...
# %         'xminortick','on',...
# %         'yminortick','on') 
# % pause


# end 

# end
# function [center] = Find_Higest_Pixel(data,yStar,xStar,d)
# % This function takes the a data cube and finds the center for a reference star
# % in each frame. The function takes the raw cube (each frame is a fits
# % file), the position of the reference star and a distance to the square.

# center_x = [];
# center_y = [];
# center   = [];
# for i = 1:length(data(1,1,:));
    
# square          = data(yStar-d:yStar+d,xStar-d:xStar+d,i); % Udsnit i hvert billede

# for k = 1:25                                % Antallet af klare pixels der skal bruges

# [pixel_i,x_max_i] = max(max(square));       % Rækken og værdien for max pixel
# [~,y_max_i]     = max(square(:,x_max_i));   % Søjlen for max pixel                   

# y_max(k)        = y_max_i;                  % Søjlerne samles i vektor  
# x_max(k)        = x_max_i;                  % Rækkerne samles i vektor
# pixel(k)        = pixel_i;                  % Max pixels samles i vektor

# % square(y_max_i,x_max_i) = 0;                % Sættes til 0 for videre beregninger

# end

# S           = sum(pixel);                   % Totale flux fra klareste pixels af udsnittet
# center_x_i  = round(1/S*dot(pixel,x_max));  % Flux-midtpunktet findes i x
# center_y_i  = round(1/S*dot(pixel,y_max));  % Flux-midtpunktet findes i y
# center_x    = [center_x center_x_i];
# center_y    = [center_y center_y_i];
# center      = [center_y' center_x'];        % center opstilles i 2 søjler


# % figure(1)
# % imshow(square,[0,1000])
# % hold on
# % plot(center_y_i,center_x_i,'b*')
# % pause
    
# end 

# end
# function [data_shift] = Shift_Frames(data,center)
# % This function shift frames. It takes a raw data cube and a 2xN vector for the center of
# % a refenrece star for the N frames.  

# offset_y    = center(:,1)-center(8,1)*ones(length(center),1);  % Offset i forhold til 
# offset_x    = center(:,2)-center(8,2)*ones(length(center),1);  % stjerne 8's koordinater

# data_shift  = [];
# for i = 1:length(center)

# % circshift forskyder billederne (yderste pixels indsættes på modsatte side). Læg mærke
# % til de to minus tegn, disse skyldes at shiftet skal være omvendt.
# data_shift_i        = circshift(data(:,:,i),[-offset_y(i) -offset_x(i)]); 
# data_shift(:,:,i)   = data_shift_i;

# end

# end
# function [i1] = Quality_Aperture(data,yStar,xStar,center,d)
# % This function calculates the quality of frames. It takes the shifted data cube,
# % the position of the reference star and a distance to the square. It uses Aperture 
# % Photometry to determined the quality of a frame.

# R_min       = round(d/5);
# R_mid       = round(d/2);
# R_max       = round(d-5);

# Quality1 = [];
# S        = [];
# B        = [];
# for i = 1:length(center)

# [yy,xx]     = meshgrid(yStar-d:yStar+d,xStar-d:xStar+d);  
# square      = data(yStar-d:yStar+d,xStar-d:xStar+d,i);      % Udsnit i hvert billede

# x           = center(i,1)+xStar-d;                  % x koordinat for stjerner
# y           = center(i,2)+yStar-d;                  % y koordinat for stjerner

# r           = sqrt((xx-x).^2 + (yy-y).^2);          % Cirklens-ligning, Radius 
                    
# Baggrund    = square(r>R_mid & r<R_max);            % Grid af flux fra baggrunden 
# B_i         = 3*median(Baggrund)-2*mean(Baggrund);  % Fluxen for baggrunden
# B           = [B B_i];

# M_i         = sum(square(r>R_min & r<R_mid)) - B_i; % Grid af flux fra ydre del af stjernen

# star        = square(r<R_min) - B_i;                % Grid af flux for indre del 
# S_i         = sum(star);                            % Fluxen fra indre del
# S           = [S S_i];

# Quality_i   = M_i/S_i;
# Quality1    = [Quality1 Quality_i];


# % % Plot af cirklerne og stjernen:
# % t           = linspace(0,2*pi,50000);
# % x1          = center(i,1);
# % y1          = center(i,2);
# % figure(1)
# % imshow(square,[0,10000])
# % hold on
# % plot(R_min*cos(t)+x1,R_min*sin(t)+y1,'g-')
# % plot(R_mid*cos(t)+x1,R_mid*sin(t)+y1,'b-')
# % plot(R_max*cos(t)+x1,R_max*sin(t)+y1,'r-')
# % plot(x1,y1,'b*')
# % 
# % pause

# end

# [~,i1]       = sort(Quality1);                                % Sortering af billeder
# % best            = round(length(data_shift(1,1,:))/100*procent); % Antallet af billeder
# % i_best          = i(1:best);                                    % Index for bedste billeder

# end
# function [i2] =Quality_Gauss_Fit(data,yStar,xStar,center,d)
# % This function calculates the quality of a frame by fitting a gauss
# % function to a reference star. The c value in the gauss function is then
# % used to determined the uality - the lower c is the sharper and well
# % defined is the peak of the star. 

# x_cut           = 1:2*d;                         
# y_cut           = 1:2*d;

# ft              = fittype('gauss1');
# opts            = fitoptions(ft);
# opts.Display    = 'Off';

# Quality2 = [];
# for i = 1:length(center)

# x               = center(i,1);                              % x koordinat for stjerner
# y               = center(i,2);                              % y koordinat for stjerner

# square          = data(yStar-d:yStar+d,xStar-d:xStar+d,i);  % Udsnit i reference billede  
    
# flux_cut_x      = square(y,x_cut);                          % flux vektor for x
# flux_cut_y      = square(y_cut,x);                          % flux vektor for y

# [xData,yData]   = prepareCurveData(x_cut,flux_cut_x);
# opts.StartPoint = [30 51 80];
# [fitresult_x]   = fit(xData,yData,ft,opts);

# [xData,yData]   = prepareCurveData(y_cut,flux_cut_y);
# opts.StartPoint = [300 340 400];
# [fitresult_y]   = fit(xData,yData,ft,opts);

# Quality_x       = fitresult_x.c1;                           % c parameteren er et udtryk 
# Quality_y       = fitresult_y.c1;                           % for bredden af gauss-kurven,
#                                                             % og denne bruges til teste                                                                                                                    
# Quality_i       = (Quality_x + Quality_y)/2;                % kvaliteten.
# Quality2        = [Quality2 Quality_i];

# end

# [~,i2]          = sort(Quality2);

# % figure(1)
# % hold on
# % plot(x_cut,flux_cut_x,'b*')
# % plot(fitresult_x,'k-')
# % set(gca,'box','on',...
# %         'xminortick','on',...
# %         'yminortick','on') 
# % xlabel('x [ pixel ]')
# % ylabel('Flux [ counts ]') 
# % grid on
# % 
# % figure(2)
# % hold on
# % plot(y_cut,flux_cut_y,'r*')
# % plot(fitresult_y,'k-')
# % set(gca,'box','on',...
# %         'xminortick','on',...
# %         'yminortick','on') 
# % xlabel('y [ pixel ]')
# % ylabel('Flux [ counts ]') 
# % grid on

# end
# function Plots(data,data_shift,i1)

# % mas          = 119.11765;

# %% Frasortering:

# bad     = [22,32,49,98,99,138,148,159,179,189,213,241,271,273,281,291,303,311,321,322,...
#            325,331,332,335,352,355,362,381,391,413,502,513,532,543,553,562,583,586,596];

# i1_bad      = ismember(i1,bad);     % Finder tal i1 som er lig med tal fra bad-vektoren. 
# i1          = [i1 i1(i1_bad)];      % Ligger disse tagerst i i1-vektoren.
# i1(i1_bad)  = [];                   % Fjerner de dobbelt eksisterende tal.
            
# % bad2     = [22,32,49,98,99,138,148,159,179,189,213,241,271,273,281,291,303,311,321,322,...
# %            325,331,332,335,352,355,362,381,391,413,502,513,532,543,553,562,583,586,596];
# % 
# % i2_bad      = ismember(i2,bad2);
# % i2          = [i2 i2(i2_bad)];
# % i2(i2_bad)  = [];

# % bad3     = [22,32,49,98,99,138,148,159,179,189,213,241,271,273,281,291,303,311,321,322,...
# %            325,331,332,335,352,355,362,381,391,413,502,513,532,543,553,562,583,586,596];
# % 
# % i3_bad      = ismember(i3,bad3);
# % i3          = [i3 i3(i3_bad)];
# % i3(i3_bad)  = [];

# %% STAR uden shift

# % 100% ikke-shifted billeder:
# frame100 = zeros(length(data(:,:,1)));
# for k = 1:length(i1) 
    
#     frame100  = frame100 + data(:,:,k);
    
# end

# % 90% shifted billeder:
# frame90 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.9 
    
#     frame90  = frame90 + data(:,:,i1(k));
    
# end

# % 80% shifted billeder:
# frame80 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.8 
    
#     frame80  = frame80 + data(:,:,i1(k));
    
# end

# % 70% shifted billeder:
# frame70 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.7 
    
#     frame70  = frame70 + data(:,:,i1(k));
    
# end

# % 60% shifted billeder:
# frame60 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.6 
    
#     frame60  = frame60 + data(:,:,i1(k));
    
# end

# % 50% shifted billeder:
# frame50 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.5 
    
#     frame50  = frame50 + data(:,:,i1(k));
    
# end

# % 40% billeder:
# frame40 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.4 
    
#     frame40  = frame40 + data(:,:,i1(k));    
# end

# % 30% billeder:
# frame30 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.3 
    
#     frame30  = frame30 + data(:,:,i1(k));    
# end

# % 20% billeder:
# frame20 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.2 
    
#     frame20  = frame20 + data(:,:,i1(k));    
# end

# % 10% billeder:
# frame10 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.1 
    
#     frame10  = frame10 + data(:,:,i1(k));    
# end

# % 5% billeder:
# frame5 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.05 
    
#     frame5  = frame5 + data(:,:,i1(k));    
# end

# % 1% billeder:
# frame1 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.01 
    
#     frame1  = frame1 + data(:,:,i1(k));    
# end


# %% i1 - Til Aperture:

# % 100% shifted billeder:
# image100 = zeros(length(data(:,:,1)));
# for k = 1:length(i1) 
    
#     image100  = image100 + data_shift(:,:,k);
    
# end

# % 90% shifted billeder:
# image90 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.9 
    
#     image90  = image90 + data_shift(:,:,i1(k));
    
# end

# % 80% shifted billeder:
# image80 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.8 
    
#     image80  = image80 + data_shift(:,:,i1(k));
    
# end

# % 70% shifted billeder:
# image70 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.7 
    
#     image70  = image70 + data_shift(:,:,i1(k));
    
# end

# % 60% shifted billeder:
# image60 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.6 
    
#     image60  = image60 + data_shift(:,:,i1(k));
    
# end

# % 50% shifted billeder:
# image50 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.5 
    
#     image50  = image50 + data_shift(:,:,i1(k));
    
# end

# % 40% billeder:
# image40 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.4 
    
#     image40  = image40 + data_shift(:,:,i1(k));    
# end

# % 30% billeder:
# image30 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.3 
    
#     image30  = image30 + data_shift(:,:,i1(k));    
# end

# % 20% billeder:
# image20 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.2 
    
#     image20  = image20 + data_shift(:,:,i1(k));    
# end

# % 10% billeder:
# image10 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.1 
    
#     image10  = image10 + data_shift(:,:,i1(k));    
# end

# % 5% billeder:
# image5 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.05 
    
#     image5  = image5 + data_shift(:,:,i1(k));    
# end

# % 1% billeder:
# image1 = zeros(length(data(:,:,1)));
# for k = 1:length(i1)*0.01 
    
#     image1  = image1 + data_shift(:,:,i1(k));    
# end


# %% i2 - Til plot af Gauss-fit

# % % 100% shifted billeder:
# % Image100 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2) 
# %     
# %     Image100  = Image100 + data_shift(:,:,k);
# %     
# % end
# % 
# % % 90% shifted billeder:
# % Image90 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2)*0.9 
# %     
# %     Image90  = Image90 + data_shift(:,:,k);
# %     
# % end
# % 
# % % 80% shifted billeder:
# % Image80 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2)*0.8 
# %     
# %     Image80  = Image80 + data_shift(:,:,k);
# %     
# % end
# % 
# % % 70% shifted billeder:
# % Image70 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2)*0.7 
# %     
# %     Image70  = Image70 + data_shift(:,:,k);
# %     
# % end
# % 
# % % 60% shifted billeder:
# % Image60 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2)*0.6 
# %     
# %     Image60  = Image60 + data_shift(:,:,k);
# %     
# % end
# % 
# % % 50% shifted billeder:
# % Image50 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2)*0.5 
# %     
# %     Image50  = Image50 + data_shift(:,:,k);
# %     
# % end
# % 
# % % 40% billeder:
# % Image40 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2)*0.4 
# %     
# %     Image40  = Image40 + data_shift(:,:,k);    
# % end
# % 
# % % 30% billeder:
# % Image30 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2)*0.3 
# %     
# %     Image30  = Image30 + data_shift(:,:,k);    
# % end
# % 
# % % 20% billeder:
# % Image20 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2)*0.2 
# %     
# %     Image20  = Image20 + data_shift(:,:,k);    
# % end
# % 
# % % 10% billeder:
# % Image10 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2)/10 
# %     
# %     Image10  = Image10 + data_shift(:,:,k);    
# % end


# %% i3 - Til COF og klareste pixel:

# % % 100% shifted billeder:
# % pic100 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3) 
# %     
# %     pic100  = pic100 + data_shift3(:,:,k);
# %     
# % end
# % 
# % % 90% shifted billeder:
# % pic90 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.9 
# %     
# %     pic90  = pic90 + data_shift3(:,:,i1(k));
# %     
# % end
# % 
# % % 80% shifted billeder:
# % pic80 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.8 
# %     
# %     pic80  = pic80 + data_shift3(:,:,i1(k));
# %     
# % end
# % 
# % % 70% shifted billeder:
# % pic70 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.7 
# %     
# %     pic70  = pic70 + data_shift3(:,:,i1(k));
# %     
# % end
# % 
# % % 60% shifted billeder:
# % pic60 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.6 
# %     
# %     pic60  = pic60 + data_shift3(:,:,i1(k));
# %     
# % end
# % 
# % % 50% shifted billeder:
# % pic50 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.5 
# %     
# %     pic50  = pic50 + data_shift3(:,:,i1(k));
# %     
# % end
# % 
# % % 40% billeder:
# % pic40 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.4 
# %     
# %     pic40  = pic40 + data_shift3(:,:,i1(k));    
# % end
# % 
# % % 30% billeder:
# % pic30 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.3 
# %     
# %     pic30  = pic30 + data_shift3(:,:,i1(k));    
# % end
# % 
# % % 20% billeder:
# % pic20 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.2 
# %     
# %     pic20  = pic20 + data_shift3(:,:,i1(k));    
# % end
# % 
# % % 10% billeder:
# % pic10 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.1 
# %     
# %     pic10  = pic10 + data_shift3(:,:,i1(k));    
# % end
# % 
# % % 5% billeder:
# % pic5 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.05 
# %     
# %     pic5  = pic5 + data_shift3(:,:,i1(k));    
# % end
# % 
# % % 1% billeder:
# % pic1 = zeros(length(data(:,:,1)));
# % for k = 1:length(i3)*0.01 
# %     
# %     pic1  = pic1 + data_shift3(:,:,i1(k));    
# % end



# %% Plot af sammenligning mellem Gauss og Aperture:
# % 
# % clf
# % 
# % x               = 480;                                     
# % y               = 290;                                      
# % x_cut           = x-10:x+10;
# % y_bar           = linspace(0,8*10^6,10^3);
# %     
# % flux_cut_100    = image100(y,x_cut);                        
# % flux_cut_90     = image90(y,x_cut);                                                 
# % flux_cut_70     = image70(y,x_cut);                                               
# % flux_cut_50     = image50(y,x_cut);                                                 
# % flux_cut_30     = image30(y,x_cut);                                                 
# % flux_cut_10     = image10(y,x_cut);                          
# % 
# % Flux_cut_100    = Image100(y,x_cut);                        
# % Flux_cut_90     = Image90(y,x_cut);                                                
# % Flux_cut_70     = Image70(y,x_cut);                                                  
# % Flux_cut_50     = Image50(y,x_cut);                                                  
# % Flux_cut_30     = Image30(y,x_cut);                                                  
# % Flux_cut_10     = Image10(y,x_cut); 
# % 
# % figure(1)
# % hold on
# % 
# % plot(x*mas,y_bar,'k-')
# % 
# % h(1) = plot(x_cut(11:end)*mas,flux_cut_100(11:end),'k-');
# % h(2) = plot(x_cut(11:end)*mas,flux_cut_90(11:end),'b-');
# % h(3) = plot(x_cut(11:end)*mas,flux_cut_70(11:end),'g-');
# % h(4) = plot(x_cut(11:end)*mas,flux_cut_50(11:end),'r-');
# % h(5) = plot(x_cut(11:end)*mas,flux_cut_30(11:end),'m-');
# % h(6) = plot(x_cut(11:end)*mas,flux_cut_10(11:end),'k--');
# % 
# % plot(x_cut(1:11)*mas,Flux_cut_100(1:11),'k-')
# % plot(x_cut(1:11)*mas,Flux_cut_90(1:11),'b-')
# % plot(x_cut(1:11)*mas,Flux_cut_70(1:11),'g-')
# % plot(x_cut(1:11)*mas,Flux_cut_50(1:11),'r-')
# % plot(x_cut(1:11)*mas,Flux_cut_30(1:11),'m-')
# % plot(x_cut(1:11)*mas,Flux_cut_10(1:11),'k--')
# % 
# % 
# % legend(h(1:6),'100 %','90 %','70 %','50 %','30 %','10 %')
# % set(gca,'FontSize',15,'box','on','xminortick','on','yminortick','on') 
# % xlabel('x (mas)')
# % ylabel('Flux (counts)') 
# % axis([5.63e4 5.817e4 0 8*10^6])
# % grid on
# % 
# % print(gcf,'/home/nicholas/Dropbox/Projekt - Obs/LaTex/Figur2.eps','-depsc')


# %% Plot af sammenligning mellem GOF og klareste pixel:

# % clf
# % 
# % x               = 480;                                     
# % y               = 290;                                      
# % x_cut           = x-10:x+10;
# % y_bar           = linspace(0,8*10^6,10^3);
# %     
# % flux_cut_100    = image100(y,x_cut);                        
# % flux_cut_90     = image90(y,x_cut);                                                 
# % flux_cut_70     = image70(y,x_cut);                                               
# % flux_cut_50     = image50(y,x_cut);                                                 
# % flux_cut_30     = image30(y,x_cut);                                                 
# % flux_cut_10     = image10(y,x_cut);                          
# % 
# % Flux_cut_100    = pic100(y,x_cut);                        
# % Flux_cut_90     = pic90(y,x_cut);                                                
# % Flux_cut_70     = pic70(y,x_cut);                                                  
# % Flux_cut_50     = pic50(y,x_cut);                                                  
# % Flux_cut_30     = pic30(y,x_cut);                                                  
# % Flux_cut_10     = pic10(y,x_cut); 
# % 
# % figure(1)
# % hold on
# % 
# % plot(x*mas,y_bar,'k-')
# % 
# % h(1) = plot(x_cut(11:end)*mas,flux_cut_100(11:end),'k-');
# % h(2) = plot(x_cut(11:end)*mas,flux_cut_90(11:end),'b-');
# % h(3) = plot(x_cut(11:end)*mas,flux_cut_70(11:end),'g-');
# % h(4) = plot(x_cut(11:end)*mas,flux_cut_50(11:end),'r-');
# % h(5) = plot(x_cut(11:end)*mas,flux_cut_30(11:end),'m-');
# % h(6) = plot(x_cut(11:end)*mas,flux_cut_10(11:end),'k--');
# % 
# % plot(x_cut(1:11)*mas,Flux_cut_100(1:11),'k-')
# % plot(x_cut(1:11)*mas,Flux_cut_90(1:11),'b-')
# % plot(x_cut(1:11)*mas,Flux_cut_70(1:11),'g-')
# % plot(x_cut(1:11)*mas,Flux_cut_50(1:11),'r-')
# % plot(x_cut(1:11)*mas,Flux_cut_30(1:11),'m-')
# % plot(x_cut(1:11)*mas,Flux_cut_10(1:11),'k--')
# % 
# % legend(h(1:6),'100 %','90 %','70 %','50 %','30 %','10 %')
# % set(gca,'FontSize',15,'box','on','xminortick','on','yminortick','on') 
# % xlabel('x (mas)')
# % ylabel('Flux (counts)') 
# % axis([5.63e4 5.817e4 0 8*10^6])
# % grid on
# % 
# % print(gcf,'/home/nicholas/Dropbox/Projekt - Obs/LaTex/Figur9.eps','-depsc')


# %% Plot til metoden: 

# % figure(1)
# % 
# % ha  = tight_subplot(2,3,[0.03 .03],[0.1 0.001],[0.01 0.01]); 
# % 
# % axes(ha(1)), imshow(data_shift(:,:,i1(end)),[0,10^(3.9)]),  xlabel('(1)'), set(gca,'Fontsize',7)
# % axes(ha(4)), imshow(data_shift(:,:,i1(1)),[0,10^(3.9)]),    xlabel('(2)'), set(gca,'Fontsize',7)
# % axes(ha(2)), imshow(frame100*(1/600),[0,10^(3.9)]),         xlabel('(3)'), set(gca,'Fontsize',7)
# % axes(ha(5)), imshow(image100*(1/600),[0,10^(3.9)]),         xlabel('(4)'), set(gca,'Fontsize',7)
# % axes(ha(3)), imshow(image50/300,[0,10^(3.9)]),              xlabel('(5)'), set(gca,'Fontsize',7) 
# % axes(ha(6)), imshow(image5/30,[0,10^(3.9)]),                xlabel('(6)'), set(gca,'Fontsize',7) 
# % 
# % axis(ha(1:6),[230 275 220 260])
# %  
# % print(gcf,'/home/nicholas/Dropbox/Projekt - Obs/LaTex/Figur1.eps','-depsc')


# %% FWHM og grafer af forskellige stjerner:

# % x               = 480;                                     
# % y               = 290;
# % x_cut           = x-20:x+20;
# % 
# % x1              = 414;                                     
# % y1              = 44;
# % x1_cut          = x1-20:x1+20;
# % 
# % x2              = 51;                                     
# % y2              = 339;
# % x2_cut          = x2-20:x2+20;
# % 
# % x3              = 426;                                     
# % y3              = 393;
# % x3_cut          = x3-20:x3+20;


# % xx               = 480;                                     
# % yy               = 293;
# % xx_cut           = xx-20:xx+20;
# % 
# % xx1              = 414;                                     
# % yy1              = 47;
# % xx1_cut          = xx1-20:xx1+20;
# % 
# % xx2              = 51;                                     
# % yy2              = 342;
# % xx2_cut          = xx2-20:xx2+20;
# % 
# % xx3              = 426;                                     
# % yy3              = 396;
# % xx3_cut          = xx3-20:xx3+20;
 
# % AperturePhotometry(image5,y1,x1,15,20)
 
 
# % 
# % flux_cut_100    = image100(y,x_cut);                        
# % flux_cut_90     = image90(y,x_cut);                        
# % flux_cut_80     = image80(y,x_cut);                         
# % flux_cut_70     = image70(y,x_cut);                         
# % flux_cut_60     = image60(y,x_cut);                         
# % flux_cut_50     = image50(y,x_cut);                         
# % flux_cut_40     = image40(y,x_cut);                         
# % flux_cut_30     = image30(y,x_cut);                         
# % flux_cut_20     = image20(y,x_cut);                         
# % flux_cut_10     = image10(y,x_cut); 
# % flux_cut_5      = image5(y,x_cut);
# % flux_cut_1      = image1(y,x_cut);
# % 
# % Flux_cut_100    = frame100(yy,xx_cut);                        
# % Flux_cut_90     = frame90(yy,xx_cut);                        
# % Flux_cut_80     = frame80(yy,xx_cut);                         
# % Flux_cut_70     = frame70(yy,xx_cut);                         
# % Flux_cut_60     = frame60(yy,xx_cut);                         
# % Flux_cut_50     = frame50(yy,xx_cut);                         
# % Flux_cut_40     = frame40(yy,xx_cut);                         
# % Flux_cut_30     = frame30(yy,xx_cut);                         
# % Flux_cut_20     = frame20(yy,xx_cut);                         
# % Flux_cut_10     = frame10(yy,xx_cut); 
# % Flux_cut_5      = frame5(yy,xx_cut);
# % Flux_cut_1      = frame1(yy,xx_cut);
# % 
# % 
# % flux1_cut_100   = image100(y1,x1_cut);                        
# % flux1_cut_90    = image90(y1,x1_cut);                        
# % flux1_cut_80    = image80(y1,x1_cut);                         
# % flux1_cut_70    = image70(y1,x1_cut);                         
# % flux1_cut_60    = image60(y1,x1_cut);                         
# % flux1_cut_50    = image50(y1,x1_cut);                         
# % flux1_cut_40    = image40(y1,x1_cut);                         
# % flux1_cut_30    = image30(y1,x1_cut);                         
# % flux1_cut_20    = image20(y1,x1_cut);                         
# % flux1_cut_10    = image10(y1,x1_cut); 
# % flux1_cut_5     = image5(y1,x1_cut);
# % flux1_cut_1     = image1(y1,x1_cut);
# % 
# % Flux1_cut_100    = frame100(yy1,xx1_cut);                        
# % Flux1_cut_90     = frame90(yy1,xx1_cut);                        
# % Flux1_cut_80     = frame80(yy1,xx1_cut);                         
# % Flux1_cut_70     = frame70(yy1,xx1_cut);                         
# % Flux1_cut_60     = frame60(yy1,xx1_cut);                         
# % Flux1_cut_50     = frame50(yy1,xx1_cut);                         
# % Flux1_cut_40     = frame40(yy1,xx1_cut);                         
# % Flux1_cut_30     = frame30(yy1,xx1_cut);                         
# % Flux1_cut_20     = frame20(yy1,xx1_cut);                         
# % Flux1_cut_10     = frame10(yy1,xx1_cut); 
# % Flux1_cut_5      = frame5(yy1,xx1_cut);
# % Flux1_cut_1      = frame1(yy1,xx1_cut);
# % 
# % 
# % flux2_cut_100    = image100(y2,x2_cut);                        
# % flux2_cut_90     = image90(y2,x2_cut);                        
# % flux2_cut_80     = image80(y2,x2_cut);                         
# % flux2_cut_70     = image70(y2,x2_cut);                         
# % flux2_cut_60     = image60(y2,x2_cut);                         
# % flux2_cut_50     = image50(y2,x2_cut);                         
# % flux2_cut_40     = image40(y2,x2_cut);                         
# % flux2_cut_30     = image30(y2,x2_cut);                         
# % flux2_cut_20     = image20(y2,x2_cut);                         
# % flux2_cut_10     = image10(y2,x2_cut); 
# % flux2_cut_5      = image5(y2,x2_cut);
# % flux2_cut_1      = image1(y2,x2_cut);
# % 
# % Flux2_cut_100    = frame100(yy2,xx2_cut);                        
# % Flux2_cut_90     = frame90(yy2,xx2_cut);                        
# % Flux2_cut_80     = frame80(yy2,xx2_cut);                         
# % Flux2_cut_70     = frame70(yy2,xx2_cut);                         
# % Flux2_cut_60     = frame60(yy2,xx2_cut);                         
# % Flux2_cut_50     = frame50(yy2,xx2_cut);                         
# % Flux2_cut_40     = frame40(yy2,xx2_cut);                         
# % Flux2_cut_30     = frame30(yy2,xx2_cut);                         
# % Flux2_cut_20     = frame20(yy2,xx2_cut);                         
# % Flux2_cut_10     = frame10(yy2,xx2_cut); 
# % Flux2_cut_5      = frame5(yy2,xx2_cut);
# % Flux2_cut_1      = frame1(yy2,xx2_cut);
# % 
# % 
# % flux3_cut_100    = image100(y3,x3_cut);                        
# % flux3_cut_90     = image90(y3,x3_cut);                        
# % flux3_cut_80     = image80(y3,x3_cut);                         
# % flux3_cut_70     = image70(y3,x3_cut);                         
# % flux3_cut_60     = image60(y3,x3_cut);                         
# % flux3_cut_50     = image50(y3,x3_cut);                         
# % flux3_cut_40     = image40(y3,x3_cut);                         
# % flux3_cut_30     = image30(y3,x3_cut);                         
# % flux3_cut_20     = image20(y3,x3_cut);                         
# % flux3_cut_10     = image10(y3,x3_cut); 
# % flux3_cut_5      = image5(y3,x3_cut);
# % flux3_cut_1      = image1(y3,x3_cut);
# % 
# % Flux3_cut_100    = frame100(yy3,xx3_cut);                        
# % Flux3_cut_90     = frame90(yy3,xx3_cut);                        
# % Flux3_cut_80     = frame80(yy3,xx3_cut);                         
# % Flux3_cut_70     = frame70(yy3,xx3_cut);                         
# % Flux3_cut_60     = frame60(yy3,xx3_cut);                         
# % Flux3_cut_50     = frame50(yy3,xx3_cut);                         
# % Flux3_cut_40     = frame40(yy3,xx3_cut);                         
# % Flux3_cut_30     = frame30(yy3,xx3_cut);                         
# % Flux3_cut_20     = frame20(yy3,xx3_cut);                         
# % Flux3_cut_10     = frame10(yy3,xx3_cut); 
# % Flux3_cut_5      = frame5(yy3,xx3_cut);
# % Flux3_cut_1      = frame1(yy3,xx3_cut);
# % 
# % xp              = [1 5 10 20 30 40 50 60 70 80 90 100];
# % 
# % ref_star        = [fwhm_new(x_cut,flux_cut_1)  fwhm_new(x_cut,flux_cut_5)...
# %                    fwhm_new(x_cut,flux_cut_10) fwhm_new(x_cut,flux_cut_20)...
# %                    fwhm_new(x_cut,flux_cut_30) fwhm_new(x_cut,flux_cut_40)...
# %                    fwhm_new(x_cut,flux_cut_50) fwhm_new(x_cut,flux_cut_60)...
# %                    fwhm_new(x_cut,flux_cut_70) fwhm_new(x_cut,flux_cut_80)...
# %                    fwhm_new(x_cut,flux_cut_90) fwhm_new(x_cut,flux_cut_100)].*mas;
# % 
# % Ref_star        = [fwhm_new(x_cut,Flux_cut_1)  fwhm_new(x_cut,Flux_cut_5)...
# %                    fwhm_new(x_cut,Flux_cut_10) fwhm_new(x_cut,Flux_cut_20)...
# %                    fwhm_new(x_cut,Flux_cut_30) fwhm_new(x_cut,Flux_cut_40)...
# %                    fwhm_new(x_cut,Flux_cut_50) fwhm_new(x_cut,Flux_cut_60)...
# %                    fwhm_new(x_cut,Flux_cut_70) fwhm_new(x_cut,Flux_cut_80)...
# %                    fwhm_new(x_cut,Flux_cut_90) fwhm_new(x_cut,Flux_cut_100)].*mas;  
# %                
# %                
# % ref1_star       = [fwhm_new(x1_cut,flux1_cut_1)  fwhm_new(x1_cut,flux1_cut_5)...
# %                    fwhm_new(x1_cut,flux1_cut_10) fwhm_new(x1_cut,flux1_cut_20)...
# %                    fwhm_new(x1_cut,flux1_cut_30) fwhm_new(x1_cut,flux1_cut_40)...
# %                    fwhm_new(x1_cut,flux1_cut_50) fwhm_new(x1_cut,flux1_cut_60)...
# %                    fwhm_new(x1_cut,flux1_cut_70) fwhm_new(x1_cut,flux1_cut_80)...
# %                    fwhm_new(x1_cut,flux1_cut_90) fwhm_new(x1_cut,flux1_cut_100)].*mas;
# % 
# % Ref1_star       = [fwhm_new(x1_cut,Flux1_cut_1)  fwhm_new(x1_cut,Flux1_cut_5)...
# %                    fwhm_new(x1_cut,Flux1_cut_10) fwhm_new(x1_cut,Flux1_cut_20)...
# %                    fwhm_new(x1_cut,Flux1_cut_30) fwhm_new(x1_cut,Flux1_cut_40)...
# %                    fwhm_new(x1_cut,Flux1_cut_50) fwhm_new(x1_cut,Flux1_cut_60)...
# %                    fwhm_new(x1_cut,Flux1_cut_70) fwhm_new(x1_cut,Flux1_cut_80)...
# %                    fwhm_new(x1_cut,Flux1_cut_90) fwhm_new(x1_cut,Flux1_cut_100)].*mas; 
# %                
# % 
# % ref2_star       = [fwhm_new(x2_cut,flux2_cut_1)  fwhm_new(x2_cut,flux2_cut_5)...
# %                    fwhm_new(x2_cut,flux2_cut_10) fwhm_new(x2_cut,flux2_cut_20)...
# %                    fwhm_new(x2_cut,flux2_cut_30) fwhm_new(x2_cut,flux2_cut_40)...
# %                    fwhm_new(x2_cut,flux2_cut_50) fwhm_new(x2_cut,flux2_cut_60)...
# %                    fwhm_new(x2_cut,flux2_cut_70) fwhm_new(x2_cut,flux2_cut_80)...
# %                    fwhm_new(x2_cut,flux2_cut_90) fwhm_new(x2_cut,flux2_cut_100)].*mas;
# % 
# % Ref2_star       = [fwhm_new(x2_cut,Flux2_cut_1)  fwhm_new(x2_cut,Flux2_cut_5)...
# %                    fwhm_new(x2_cut,Flux2_cut_10) fwhm_new(x2_cut,Flux2_cut_20)...
# %                    fwhm_new(x2_cut,Flux2_cut_30) fwhm_new(x2_cut,Flux2_cut_40)...
# %                    fwhm_new(x2_cut,Flux2_cut_50) fwhm_new(x2_cut,Flux2_cut_60)...
# %                    fwhm_new(x2_cut,Flux2_cut_70) fwhm_new(x2_cut,Flux2_cut_80)...
# %                    fwhm_new(x2_cut,Flux2_cut_90) fwhm_new(x2_cut,Flux2_cut_100)].*mas; 
# %      
# %                
# % ref3_star       = [fwhm_new(x3_cut,flux3_cut_1)  fwhm_new(x3_cut,flux3_cut_5)...
# %                    fwhm_new(x3_cut,flux3_cut_10) fwhm_new(x3_cut,flux3_cut_20)...
# %                    fwhm_new(x3_cut,flux3_cut_30) fwhm_new(x3_cut,flux3_cut_40)...
# %                    fwhm_new(x3_cut,flux3_cut_50) fwhm_new(x3_cut,flux3_cut_60)...
# %                    fwhm_new(x3_cut,flux3_cut_70) fwhm_new(x3_cut,flux3_cut_80)...
# %                    fwhm_new(x3_cut,flux3_cut_90) fwhm_new(x3_cut,flux3_cut_100)].*mas;
# % 
# % Ref3_star       = [fwhm_new(x3_cut,Flux3_cut_1)  fwhm_new(x3_cut,Flux3_cut_5)...
# %                    fwhm_new(x3_cut,Flux3_cut_10) fwhm_new(x3_cut,Flux3_cut_20)...
# %                    fwhm_new(x3_cut,Flux3_cut_30) fwhm_new(x3_cut,Flux3_cut_40)...
# %                    fwhm_new(x3_cut,Flux3_cut_50) fwhm_new(x3_cut,Flux3_cut_60)...
# %                    fwhm_new(x3_cut,Flux3_cut_70) fwhm_new(x3_cut,Flux3_cut_80)...
# %                    fwhm_new(x3_cut,Flux3_cut_90) fwhm_new(x3_cut,Flux3_cut_100)].*mas; 
# %  
# %                
# %                
# % figure(1)
# % 
# % subplot(1,3,[1 2])
# % hold on 
# % plot(xp,ref_star,'--ro','LineWidth',1,'MarkerEdgeColor','k',...
# %                                     'MarkerFaceColor','r',...
# %                                     'MarkerSize',6)
# %                                                                
# % plot(xp,ref1_star,'--mo','LineWidth',1,'MarkerEdgeColor','k',...
# %                                     'MarkerFaceColor','m',...
# %                                     'MarkerSize',6)
# %                                  
# % plot(xp,ref2_star,'--go','LineWidth',1,'MarkerEdgeColor','k',...
# %                                     'MarkerFaceColor','g',...
# %                                     'MarkerSize',6) 
# % 
# % plot(xp,ref3_star,'--bo','LineWidth',1,'MarkerEdgeColor','k',...
# %                                     'MarkerFaceColor','b',...
# %                                     'MarkerSize',6)                                
# % 
# % plot(xp,Ref_star,'--rs','LineWidth',1,'MarkerEdgeColor','k',...
# %                                     'MarkerFaceColor','r',...
# %                                     'MarkerSize',6) 
# %                                 
# % plot(xp,Ref1_star,'--ms','LineWidth',1,'MarkerEdgeColor','k',...
# %                                     'MarkerFaceColor','m',...
# %                                     'MarkerSize',6)                                
# %                                
# % plot(xp,Ref2_star,'--gs','LineWidth',1,'MarkerEdgeColor','k',...
# %                                     'MarkerFaceColor','g',...
# %                                     'MarkerSize',6)
# %                                 
# % plot(xp,Ref3_star,'--bs','LineWidth',1,'MarkerEdgeColor','k',...
# %                                     'MarkerFaceColor','b',...
# %                                     'MarkerSize',6)    
# %                                 
# % set(gca,'Fontsize',12,'box','on','xminortick','on','yminortick','on') 
# % xlabel('Procent')
# % ylabel('FWHM (mas)') 
# % % legend('With shift','Without shift','Location','NorthWest')
# % axis([0 102 275 825])
# % grid on
# % 
# % print(gcf,'/home/nicholas/Dropbox/Projekt - Obs/LaTex/Figur3.eps','-depsc')


# %% Fejl i billederne:

# % figure(1)
# % 
# % subplot(1,3,1)
# % imshow(data(:,:,88),[0,6500])
# % set(gca,'Fontsize',7)
# % axis([83 103 215 235])
# % 
# % subplot(1,3,2)
# % imshow(data(:,:,241),[0,6500])
# % set(gca,'Fontsize',7)
# % axis([299 319 370 390])
# % 
# % subplot(1,3,3)
# % imshow(data(:,:,332),[0,6500])
# % % colorbar
# % set(gca,'Fontsize',7)
# % axis([208 228 298 318])
# % 
# % print(gcf,'/home/nicholas/Dropbox/Projekt - Obs/LaTex/Figur5.eps','-depsc')
# % 
# % figure(2)
# % 
# % ax(1) = subplot(1,3,1)
# % imshow(data(:,:,87),[0,6500])
# % set(gca,'Fontsize',7)
# % axis([83 103 215 235])
# % 
# % ax(2) = subplot(1,3,2)
# % imshow(data(:,:,240),[0,6500])
# % set(gca,'Fontsize',7)
# % axis([299 319 370 390])
# % 
# % ax(3) = subplot(1,3,3)
# % imshow(data(:,:,331),[0,6500])
# % set(gca,'Fontsize',7)
# % axis([208 228 298 318])
# % 
# % print(gcf,'/home/nicholas/Dropbox/Projekt - Obs/LaTex/Figur6.eps','-depsc')

# %% Plot til Aperture metoden: 

# % % Plot af cirklerne og stjernen:
# % t           = linspace(0,2*pi,50000);
# % x1          = 480;
# % y1          = 290;
# % R_min       = round(20/5);
# % R_mid       = round(20/2);
# % R_max       = round(20-5);
# % 
# % figure(1)
# % imshow(data_shift(:,:,i1(1)),[0,6000])
# % hold on
# % plot(R_min*cos(t)+x1,R_min*sin(t)+y1,'g-','LineWidth',2)
# % plot(R_mid*cos(t)+x1,R_mid*sin(t)+y1,'b-','LineWidth',2)
# % plot(R_max*cos(t)+x1,R_max*sin(t)+y1,'r-','LineWidth',2)
# % axis([x1-20 x1+20 y1-20 y1+20])
# % 
# % print(gcf,'/home/nicholas/Dropbox/Projekt - Obs/LaTex/Figur7.eps','-depsc')


# %% 3 billeder af hele feltet:

# % figure(1)
# % 
# % ha  = tight_subplot(1,3,[0.03 .03],[0.1 0.2],[0.01 0.01]); 
# % 
# % axes(ha(1)), imshow(frame100 /600,[0,10^(3.9)])
# % axes(ha(2)), imshow(image100 /600,[0,10^(3.9)])
# % axes(ha(3)), imshow(image5 /30,[0,10^(3.9)])
# % 
# % axis(ha(1:3),[100 430 90 400])
# % 
# % print(gcf,'/home/nicholas/Dropbox/Projekt - Obs/LaTex/Figur8.eps','-depsc')


# %% Signal-to-Noise plot:


# % SN      = [354 756.3 1056 1484 1816 2074 2309 2506 2703 2864 3016 3135];
# % 
# % x3              = 426;                                     
# % y3              = 393;
# % x3_cut          = x3-20:x3+20;
# % 
# % AperturePhotometry(data_shift(:,:,1),y3,x3,7,13)

# % 
# % flux3_cut_100    = image100(y3,x3_cut);                        
# % flux3_cut_90     = image90(y3,x3_cut);                        
# % flux3_cut_80     = image80(y3,x3_cut);                         
# % flux3_cut_70     = image70(y3,x3_cut);                         
# % flux3_cut_60     = image60(y3,x3_cut);                         
# % flux3_cut_50     = image50(y3,x3_cut);                         
# % flux3_cut_40     = image40(y3,x3_cut);                         
# % flux3_cut_30     = image30(y3,x3_cut);                         
# % flux3_cut_20     = image20(y3,x3_cut);                         
# % flux3_cut_10     = image10(y3,x3_cut); 
# % flux3_cut_5      = image5(y3,x3_cut);
# % flux3_cut_1      = image1(y3,x3_cut);
# % 
# % ref3_star       = [fwhm_new(x3_cut,flux3_cut_1)  fwhm_new(x3_cut,flux3_cut_5)...
# %                    fwhm_new(x3_cut,flux3_cut_10) fwhm_new(x3_cut,flux3_cut_20)...
# %                    fwhm_new(x3_cut,flux3_cut_30) fwhm_new(x3_cut,flux3_cut_40)...
# %                    fwhm_new(x3_cut,flux3_cut_50) fwhm_new(x3_cut,flux3_cut_60)...
# %                    fwhm_new(x3_cut,flux3_cut_70) fwhm_new(x3_cut,flux3_cut_80)...
# %                    fwhm_new(x3_cut,flux3_cut_90) fwhm_new(x3_cut,flux3_cut_100)].*mas;
# % 
# % 
# % xp      = [1 5 10 20 30 40 50 60 70 80 90 100];
# % 
# % figure(1)
# % hold on
# % 
# % [AX,H1,H2] = plotyy(xp,ref3_star,xp,SN);
# % 
# % ylabel(AX(1),'FWHM','FontSize',17)             % label left y-axis
# % ylabel(AX(2),'Signal-to-Noise','FontSize',17)  % label right y-axis
# % xlabel(AX(2),'Procent','FontSize',17)          % label x-axis
# % 
# % set(H1,'linestyle','--','marker','o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',7) 
# % set(H2,'linestyle','--','marker','o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',7)
# % set(AX,'FontSize',17,'box','on','xminortick','on')%,'yminortick','on')
# % 
# % grid on
# % 
# % set(AX(1),'YLim',[300 550])
# % set(AX(1),'YTick',[300:50:550])
# % set(AX(2),'YLim',[0 3500])
# % set(AX(2),'YTick',[0:500:3500])
# % 
# % print(gcf,'/home/nicholas/Dropbox/Projekt - Obs/LaTex/Figur10.eps','-depsc')


# %% 2 slutbilleder over M13:

# % x               = 480;                                     
# % y               = 290;
# % 
# % x1              = 414;                                     
# % y1              = 44;
# % 
# % x2              = 51;                                     
# % y2              = 339;
# % 
# % x3              = 426;                                     
# % y3              = 393;
# % 
# % t           = linspace(0,2*pi,50000);
# % R_min       = 10;
# % 
# % l1          = 100:430;
# % l2          = 90:400;
# % 
# % r1          = 230:275;
# % r2          = 220:260;
# % 
# % figure(1)
# % % Sidste billede
# % imshow(image5/30,[0,10^(3.9)])
# % hold on
# % 
# % plot(l1,90,'y-')
# % plot(l1,400,'y-')
# % plot(100,l2,'y-')
# % plot(430,l2,'y-')
# % 
# % plot(r1,220,'c-')
# % plot(r1,260,'c-')
# % plot(230,r2,'c-')
# % plot(275,r2,'c-')
# % 
# % plot(R_min*cos(t)+x,R_min*sin(t)+y,'r-')
# % plot(R_min*cos(t)+x1,R_min*sin(t)+y1,'m-')
# % plot(R_min*cos(t)+x2,R_min*sin(t)+y2,'g-')
# % plot(R_min*cos(t)+x3,R_min*sin(t)+y3,'b-')
# % 
# % print(gcf,'/home/nicholas/Dropbox/Projekt - Obs/LaTex/Figur11.eps','-depsc')    

# %% Værdier for FWHM:

# % 1% billeder:
# % Image1 = zeros(length(data(:,:,1)));
# % for k = 1:length(i2)*0.01 
# %     
# %     Image1  = Image1 + data_shift(:,:,k);    
# % end

# % x               = 480;                                     
# % y               = 290;                                      
# % x_cut           = x-10:x+10;
# %                                                 
# % flux_cut_1      = image1(y,x_cut);
# % flux_cut_100    = image100(y,x_cut);
# % flux_cut_100_x  = frame100(y,x_cut);
# % 
# % a   = fwhm_new(x_cut,flux_cut_1)
# % b   = fwhm_new(x_cut,flux_cut_100)
# % c   = fwhm_new(x_cut,flux_cut_100_x)
# % 
# % c/b
# % b/a

# %%%%%%%%%%%% for the very faint star:

# % x3              = 426;                                     
# % y3              = 393;
# % x3_cut          = x3-20:x3+20;
# % 
# % flux3_cut_1      = image1(y3,x3_cut);
# % flux3_cut_100    = image100(y3,x3_cut);                        
# % Flux3_cut_100    = frame100(y3,x3_cut); 
# % 
# % a   = fwhm_new(x3_cut,flux3_cut_1)
# % b   = fwhm_new(x3_cut,flux3_cut_100)
# % c   = fwhm_new(x3_cut,Flux3_cut_100)
# % 
# % d1 = c/b;
# % d2 = b/a;
# % 
# % d1+d2

# %% EFTER PROJEKTET:

# x              = 425;                                     
# y              = 393;
# x_cut          = x-20:x+20;

# % figure(1)
# % imshow(pic5/30,[0,10^4])
# % hold on
# % plot(x,y,'b*')

# % figure(1)
# % ha  = tight_subplot(2,2,[0.03 .03],[0.1 0.001],[0.01 0.01]); 
# % axes(ha(1)), imshow(frame100/600,[0,10^(3.9)])
# % axes(ha(2)), imshow(frame100/600,[0,10^(3.9)])   
# % axes(ha(3)), imshow(image100/600,[0,10^(3.9)])         
# % axes(ha(4)), imshow(pic100/600,[0,10^(3.9)])         
# % axis(ha(1:4),[230 275 220 260])

# Flux_cut_100    = frame100(y,x_cut); 
# flux_cut_100    = image100(y,x_cut);
# flux_cut_1      = image1(y,x_cut);

# a1   = fwhm_new(x_cut,flux_cut_1)
# b1   = fwhm_new(x_cut,flux_cut_100)
# c1   = fwhm_new(x_cut,Flux_cut_100)

# d_1 = c1/b1
# d_2 = b1/a1
# d_3 = c1/a1

# % flux_cut_pix_100    = pic100(y,x_cut);
# % flux_cut_pix_1      = pic1(y,x_cut);

# % a2   = fwhm_new(x_cut,flux_cut_pix_1)
# % b2   = fwhm_new(x_cut,flux_cut_pix_100);
# % c2   = fwhm_new(x_cut,Flux_cut_100);

# % d_2 = c2/b2
# % d_2 = b2/a2
# % d_2 = c2/a2


# end






