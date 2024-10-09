function retv=Fcn_Bruker_Method(pname, fname, vname)
% Fcn_Bruker_Method Read bruker configuration files: sacqp and method fiiles 
%
%   Input: 
%       pname: path to the file
%       fname: configuration file 
%       vname: varialbe name
% 
%   output: 
%       retv: parameter value

if pname(end)~='/'
    pname=[pname '/']; 
end

fid=fopen([pname fname]);
retv=-999999999; % nothing found

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    sp0=strfind(tline, ['##$' vname]);
    if length(sp0) > 0 % found
        sp=strfind(tline, '(');
        if length(sp) > 0 % array
            sp2=strfind(tline, ')');
            narray=str2num(tline(sp+1:sp2-1));
            if length(narray) == 1
                retv=zeros(1, narray);
                pt=1;
                while pt <= narray
                    tline=fgetl(fid);
                    if tline(1)=='<'
                        retv=tline(2:end-1);
                        break;
                    end
                    sp=strfind(tline, ' ');
                    sp=[1 sp];
                    if sp(end) < length(tline), sp=[sp length(tline)]; end
                    if sp(end) == 1,
                        retv(pt)=str2num(tline);
                        pt=pt+1;
                    else
                        for d=1:length(sp)-1
                            temp=str2num(tline(sp(d):sp(d+1)));
                            if length(temp)>0
                                retv(pt)=temp;
                                pt=pt+1;
                            else
                                retv=tline(sp(d):sp(d+1)); % it's a string!
                                pt=narray+1; % ignore the rest
                                break;
                            end % if
                        end % for d=1
                    end % if sp(end) == 1
                end % while
            else % if dimension > 1
                retv=zeros([narray]);
                ai=ones(1, length(narray)); % array index
                npt=prod(narray); % total number of points
                pt=1;
                while pt <= npt
                    tline=fgetl(fid);
                    if tline(1)=='<'
                        retv=tline(2:end-1);
                        break;
                    end
                    sp=strfind(tline, ' ');
                    sp=[1 sp];
                    if sp(end) < length(tline), sp=[sp length(tline)]; end
                    if sp(end) == 1,
                        ind=bruker_method_aind(ai, narray);
                        retv(ind)=str2num(tline);
                        pt=pt+1;
                        ai=bruker_method_inc_ai(ai, narray);
                        %[ind ai]
                    else
                        for d=1:length(sp)-1
                            temp=str2num(tline(sp(d):sp(d+1)));
                            if length(temp)>0
                                ind=bruker_method_aind(ai, narray);
                                retv(ind)=temp;
                                pt=pt+1;
                                ai=bruker_method_inc_ai(ai, narray);
                                %[ind ai]
                            else
                                retv=tline(sp(d):sp(d+1)); % it's a string! shouldn't be though.
                                pt=npt+1; % ignore the rest
                                break;
                            end % if
                        end % for d=1
                    end % if sp(end) == 1
                end % while
            end % dimension > 1
            
        else % if not array
            retv=str2num(tline(sp0(1)+length(['##$' vname])+1:end));
            if length(retv)==0 % not a number string
                retv=tline(sp0(1)+length(['##$' vname])+1:end);
            end
        end
        break; % completed
    end
end
fclose(fid);

return

function ai=bruker_method_inc_ai(ai, narray)
%[ai narray]
na=length(narray);
ai(na)=ai(na)+1;
for n=na:-1:2
    if ai(n) > narray(n)
        ai(n) = 1;
        ai(n-1) = ai(n-1)+1;
    end
end
if ai(1) > narray(1), ai(1)=1; end
return

function ind=bruker_method_aind(ai, narray)
na=length(narray);
ind=1;
for n=1:na
    if n==1
        ind=ind+(ai(n)-1);
    else
        ind=ind+(ai(n)-1)*prod(narray(1:n-1));
    end
end
%[ind]
return
