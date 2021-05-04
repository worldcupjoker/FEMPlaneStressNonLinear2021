% Function: Read input file for the specified string
% Input   : [string]   specifies data to be read; 
%            [filename] specifies name of .inp file to be read
% Output  : 1D or 2D array of data

function A = readinp(string,filename)

fid = fopen(filename, 'r');
A = [];
while ~feof(fid)
    s = fgetl(fid);
    if (findstr(s, string)==1)
        i=1;
        while ~feof(fid)
            s = fgetl(fid);
            if (size(findstr(s, '*'),1)~=0)
                break;
            end
            if (size(strfind(string, 'set'),1)>0)
                A=[A sscanf(s,'%f%*s')'];
            else
                A(i,:) = sscanf(s,'%f%*s')';
                i=i+1;
            end
        end
        break;
    end
end
fclose(fid);


