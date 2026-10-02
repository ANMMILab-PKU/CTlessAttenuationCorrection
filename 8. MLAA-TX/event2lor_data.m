function [sino] = event2lor(varargin)
% MODIFIED: Accepts crystal_lut array directly to avoid repeated I/O

%% Initialization
if isnumeric(varargin{1})
  format='array';
else
  filename=varargin{1};
  format=filename((length(filename)-2):length(filename));
end

% Check if the second argument is a filename (char) or the data array
lut_input = varargin{2};
if ischar(lut_input) || isstring(lut_input)
    crystal_position = readlut(lut_input);
else
    crystal_position = lut_input;
end

 if nargin == 3
   PET_geom=varargin{3};
   cry_num_tr = PET_geom(:,1);
   cry_num_ax = PET_geom(:,2);
   mod_num_tr = PET_geom(:,3);   
   rsector_num = PET_geom(:,4);
   ring_num = PET_geom(:,5);
 elseif nargin <3
    error('Not enough inputs!')
 else
    error('Too many inputs!')
 end

[layer_num,temp]=size(PET_geom);
M=cry_num_ax.*ring_num;
L=cry_num_tr.*mod_num_tr.*rsector_num;

%% sino-matrix

switch format
    case 'mat'  
       load(filename)
       data=coin_final_result_all;
       ID=zeros(length(data),3);
       % Calculation of IDs logic retained...
       ID(:,2) =  ( ( (ring_num(1)-1-data(:,6))*cry_num_ax(1)  + floor(data(:,4)/cry_num_ax(1)) ) * cry_num_tr(1) * rsector_num(1) ...
               +  mod(rsector_num(1)-data(:,5),rsector_num(1))* cry_num_tr(1) ...
               +  cry_num_tr(1)-mod(data(:,4),cry_num_ax(1))-1) ...
               .* (data(:,4) < (cry_num_tr(1)*cry_num_ax(1))); 
       ID(:,3) =  ( ( (ring_num(1)-1-data(:,12))*cry_num_ax(1)  + floor(data(:,10)/cry_num_ax(1)) ) * cry_num_tr(1) * rsector_num(1) ...
               +  mod(rsector_num(1)-data(:,11),rsector_num(1))* cry_num_tr(1) ...
               +  cry_num_tr(1)-mod(data(:,10),cry_num_ax(1))-1) ...
               .* (data(:,10) < (cry_num_tr(1)*cry_num_ax(1)));  

    case {'cdf','array'}
      switch format
          case 'cdf'
          cdffile =fopen(filename,'r','l');
          ID = fread(cdffile,'int32');
          ID = reshape(ID,3,[])';

          case 'array'
          ID = varargin{1};
          ID = [zeros(length(ID),1) ID];
      end    
    otherwise
    error('Unknown format! Unable to process!')
end 

   % Generate coordinates using the passed LUT array
   data=zeros(length(ID),12);
   
   % Ensure indices are within bounds
   data(:,1:3) = crystal_position(ID(:,2)+1,1:3);
   data(:,7:9) = crystal_position(ID(:,3)+1,1:3);

   % 3D_ID
   sino=zeros(length(ID),10);
   sino(:,1)=ID(:,2);
   sino(:,4)=ID(:,3);
     
   % ring_ID
   sino(:,2) = floor(sino(:,1)/L(1)) .* (sino(:,1)<M(1)*L(1));
   sino(:,5) = floor(sino(:,4)/L(1)) .* (sino(:,4)<M(1)*L(1));

   % 2D_ID
   sino(:,3) = mod(sino(:,1),L(1)).* (sino(:,1)<(M(1)*L(1)));
   sino(:,6) = mod(sino(:,4),L(1)).* (sino(:,4)<(M(1)*L(1)));

   % r
   r=abs( (data(:,2).*data(:,7)-data(:,8).*data(:,1))./sqrt((data(:,8)-data(:,2)).^2+(data(:,7)-data(:,1)).^2) );
   r=r.*sign( (abs(sino(:,3)-sino(:,6)) < (0.5 *L(1)))-0.5 );
   sino(:,7)=round(100*r)/100;

   % phi   
   u=[1 ; 0]; 
   cosine=abs(  [data(:,1)-data(:,7)  data(:,2)-data(:,8)]*u ...
        ./ sqrt((data(:,1)-data(:,7)).^2+(data(:,2)-data(:,8)).^2)  );
   cosine  =  cosine.*sign( (min(sino(:,3),sino(:,6)) < (0.25 * L(1) - 1)) -0.5 );
   phi=acosd(cosine);
   sino(:,8)=real(round(10*phi)/10); 

   % theta
   v=[0 ; 0 ; 1]; 
   cosine2=abs(  [data(:,1)-data(:,7)  data(:,2)-data(:,8)  data(:,3)-data(:,9)]*v ...
        ./ sqrt((data(:,1)-data(:,7)).^2+(data(:,2)-data(:,8)).^2+(data(:,3)-data(:,9)).^2)  );
   theta=90-acosd(cosine2);
   sino(:,9)=real(round(10*theta)/10);

   % delta_z
   sino(:,10)= sino(:,2) + sino(:,5) - (M(1)-1)/2 ;

end
