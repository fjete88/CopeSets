function sImg = mySmooth(Img,smo)

%__________________________________________________________________________
% Smooths a signal and normalizes to max unity
% Input:
%	Img    - Signal field (can be 2D or 3D)
%   smo    - Smoothing FWHM
%
% Output:
%   sImg   - Smoothed signal that has been renormalized to max unity
%_________________________________________________________________________

sImg = zeros(size(Img));

if any(smo)
  spm_smooth(double(Img),sImg,smo);
  sImg = sImg/max(sImg(:));
else
  sImg = Img/max(double(Img(:)));
end

Img(Img(:)<0.05) = 0;

return