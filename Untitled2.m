
fileID = fopen("D:/Download/Chrome/student",'r','n','UTF-8');
students = {}; % creat a students cell
n = 1; 
studentchar = fgetl(fileID); % read the 1st line
studentcell = strsplit(studentchar,','); % split student info from char to cell
student = struct(); % prepare a empty struct and studentcell features will be in
while isa(studentcell,'cell')
    student.name = studentcell{1};
    student.ID = studentcell{2};
    student.home = studentcell{3};
    student.score = str2double(studentcell{4});
    students(n) = {student};
    n = n+1;
    studentchar = fgetl(fileID);
    if isa(studentchar,'char') % add this iteration, because we know the end of fgetl is -1, which can not be split
        studentcell = strsplit(studentchar,',');
    else
        studentcell = studentchar;
    end
end