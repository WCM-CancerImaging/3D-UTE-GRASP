function []=Fcn_write_analyze(img, dim, pixdim, pname, fname, dtype, origin)
% Fcn_write_analyze Write data matrix as analyze file 
%
%   Input: 
%       img: image matrix
%       dim: image dimension
%       pixdim: image pixel dimention
%       pname: path name
%       fname: file name
%       dtype: data type
%       origin: origin
% 
%   Output: 
%       analyze images are saved to local directory

if ~exist('dtype'), dtype='float32'; end
switch dtype
    case 'uint8'
        datatype=2;
        bitpix=8;
    case 'int16'
        datatype=4;
        bitpix=16;
    case 'int32'
        datatype=8;
        bitpix=32;
    case 'float32'
        datatype=16;
        bitpix=32;
    case 'complex'
        datatype=32;
        bitpix=64;
    case 'float64'
        datatype=64;
        bitpix=64;
    otherwise
        datatype=16;
        dtype='float32';
        bitpix=64;
end
        
if ~exist('fname'),
    [fname,pname] = uiputfile(['*.img'],'Enter filename(without .img) to save');
    if isempty(fname), return; end
end
six=findstr(fname, '.');
sixn=length(six);
if sixn >0 
   fcore=fname(1:six(end)-1);
else
   fcore=fname;
end
%if pname(end) ~= '\', pname=[pname '\']; end
if pname(end) ~= '/', pname=[pname '/']; end % specifically for HPC linux format

img2=zeros(dim(2), dim(1), dim(3));

for i=1:dim(3)
    img2(:,:,i)=(rot90(img(:,:,i), 3));
end
fid=fopen([pname fcore '.img'], 'w');
fwrite(fid, img2, dtype);
fclose(fid);
    

% write header
fid=fopen([pname fcore '.hdr'], 'w');
fwrite(fid,348,'int32');             %hdr.hk.sizeof_hdr
fwrite(fid,zeros(10,1),'uint8');     %hdr.hk.data_type
fwrite(fid,zeros(18, 1),'uint8');    %hdr.hk.db_name
fwrite(fid,0,'int32');               %hdr.hk.extents
fwrite(fid,0,'int16');               %hdr.hk.session_error
fwrite(fid,114,'int8');              %hdr.hk.regular
fwrite(fid,0,'int8');                %hdr.hk.hkey_un0

temp=zeros(8,1);
temp(1)=length(dim);
temp(2)=dim(2);
temp(3)=dim(1);
temp(4)=dim(3);
fwrite(fid,temp,'int16');            %hdr.dime.dim
fwrite(fid,0,'int16');               %hdr.dime.unused8
fwrite(fid,0,'int16');               %hdr.dime.unused9
fwrite(fid,0,'int16');               %hdr.dime.unused10
fwrite(fid,0,'int16');               %hdr.dime.unused11
fwrite(fid,0,'int16');               %hdr.dime.unused12
fwrite(fid,0,'int16');               %hdr.dime.unused13
fwrite(fid,0,'int16');               %hdr.dime.unused14
fwrite(fid,uint16(datatype),'int16');              %hdr.dime.datatype
fwrite(fid,uint16(bitpix),'int16');              %hdr.dime.bitpix
fwrite(fid,0,'int16');               %hdr.dime.dim_un0
temp=zeros(8,1);
temp(2)=pixdim(2);
temp(3)=pixdim(1);
temp(4)=pixdim(3);
fwrite(fid,temp,'float32');          %hdr.dime.pixdim
fwrite(fid,0,'float32');             %hdr.dime.vox_offset=
fwrite(fid,0,'float32');             %hdr.dime.funused1=
fwrite(fid,0,'float32');             %hdr.dime.funused2=
fwrite(fid,0,'float32');             %hdr.dime.funused3=
fwrite(fid,0,'float32');             %hdr.dime.cal_max=
fwrite(fid,0,'float32');             %hdr.dime.cal_min=
fwrite(fid,0,'float32');             %hdr.dime.compressed=
fwrite(fid,0,'float32');             %hdr.dime.verified=
fwrite(fid,0,'int32');               %hdr.dime.glmax=
fwrite(fid,0,'int32');               %hdr.dime.glmin=

fwrite(fid,zeros(80,1),'uchar');     %hdr.hist.descrip=
fwrite(fid,zeros(24,1),'uchar');     %hdr.hist.aux_file=
fwrite(fid,0,'int8');                %hdr.hist.orient=
temp=zeros(10,1);
if ~exist('origin'),
    temp(1)=mod(dim(2)/2, 256);
    temp(2)=floor(dim(2)/2/256);
    temp(3)=mod(dim(1)/2, 256);
    temp(4)=floor(dim(1)/2/256);
    temp(5)=mod(dim(3)/2, 256);
    temp(6)=floor(dim(3)/2/256);
else
    temp(1)=mod(origin(1), 256); % x
    temp(2)=floor(origin(1)/256); % x
    temp(3)=mod(origin(2), 256); % y 
    temp(4)=floor(origin(2)/256); % y 
    temp(5)=mod(origin(3), 256); % z
    temp(6)=floor(origin(3)/256); % z
end

fwrite(fid,temp,'uint8');             %hdr.hist.originator=
temp=zeros(10,1);
fwrite(fid,temp,'int8');             %hdr.hist.generated=
fwrite(fid,temp,'int8');             %hdr.hist.scannum=
fwrite(fid,temp,'int8');             %hdr.hist.patient_id=
fwrite(fid,temp,'int8');             %hdr.hist.exp_data=
fwrite(fid,temp,'int8');             %hdr.hist.exp_time=
fwrite(fid,zeros(3,1),'int8');       %hdr.hist.hist_un0=
fwrite(fid,0,'int32');               %hdr.hist.views=
fwrite(fid,0,'int32');               %hdr.hist.vols_added=
fwrite(fid,0,'int32');               %hdr.hist.start_field=
fwrite(fid,0,'int32');               %hdr.hist.field_skip=
fwrite(fid,0,'int32');               %hdr.hist.omax=
fwrite(fid,0,'int32');               %hdr.hist.omin=
fwrite(fid,0,'int32');               %hdr.hist.smax=
fwrite(fid,0,'int32');               %hdr.hist.smin=
fclose(fid);
% end of reading header
