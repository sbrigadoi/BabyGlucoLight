function [nirs_file_number_matlab,nirs_file_number,N_seperate_recordings,rec_length] = get_nirs_file_number_matlab(subjectN)


%each subject needs a new case
%with new N_seperate_recordings and new rec_length
switch subjectN

    case 5
        N_seperate_recordings = 1;
        rec_length = [113];

        nirs_file_number_matlab_1 = [1 ;
            10;
            (100:1:109)';
            11;
            (110:1:113)';
            12;
            13;
            14;
            (15:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:79)';
            8;
            (80:1:89)';
            9;
            (90:1:99)';
            0];


    case 9
    N_seperate_recordings = 2;
        rec_length = [77 29];

        nirs_file_number_matlab_1 = [1 ;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:77)';
            8;
            9;
            0];

            nirs_file_number_matlab_2 = [1 ;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            4;
            5;
            6;
            7;
            8;
            9;
            0];



    case 10
        N_seperate_recordings = 1;
        rec_length = [91];
        nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:79)';
            8;
            (80:1:89)';
            9;
            (90:1:91)';
            0];


    case 11
        N_seperate_recordings = 1;
        rec_length = [30];
        nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            30;
            4;
            5;
            6;
            7;
            8;
            9;
            0];

    case 17
        N_seperate_recordings = 1;
        rec_length = [45];
        nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:45)';
            5;
            6;
            7;
            8;
            9;
            0];


     case 28
        N_seperate_recordings = 1;
        rec_length = [56];
        nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:56)';
            6;
            7;
            8;
            9;
            0];

    case 34
        N_seperate_recordings = 1;
        rec_length = [56];
        nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:56)';
            6;
            7;
            8;
            9;
            0];

    case 35
        N_seperate_recordings = 1;
        rec_length = [78];
        nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:78)';
            8;
            9;
            0];


    case 36
    
        N_seperate_recordings = 1;
        rec_length = [159];

        nirs_file_number_matlab_1 = [1 ;
            10;
            (100:1:109)';
            11;
            (110:1:119)';
            12;
            (120:1:129)';
            13;
            (130:1:139)';
            14;
            (140:1:149)';
            15;
            (150:1:159)';
            (16:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:79)';
            8;
            (80:1:89)';
            9;
            (90:1:99)';
            0];

    case 39
        N_seperate_recordings = 1;
        rec_length = [50];
        nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            50;
            6;
            7;
            8;
            9;
            0];

    case 41
        N_seperate_recordings = 4;
        rec_length = [18 52 18 26];

        nirs_file_number_matlab_1 = [1;
            (10:1:18)';
            2;
            3;
            4;
            5;
            6;
            7;
            8;
            9;
            0];

         nirs_file_number_matlab_2 = [1;
             (10:1:19)';
             2;
             (20:1:29)';
             3;
             (30:1:39)';
             4;
             (40:1:49)';
             5;
             (50:1:52)';
             6;    
             7;
             8;
             9;
             0];

            nirs_file_number_matlab_3 = [1;
            (10:1:18)';
            2;
            3;
            4;
            5;
            6;
            7;
            8;
            9;
            0];

            nirs_file_number_matlab_4 = [1;
             (10:1:19)';
             2;
             (20:1:26)';
             3; 
             4;
             5; 
             6;    
             7;
             8;
             9;
             0];
       

    case 43
    
        N_seperate_recordings = 1;
        rec_length = [117];

        nirs_file_number_matlab_1 = [1 ;
            10;
            (100:1:109)';
            11;
            (110:1:117)';
            12;
            13;
            14;
            15;
            (16:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:79)';
            8;
            (80:1:89)';
            9;
            (90:1:99)';
            0];

    case 44
        N_seperate_recordings = 1;
        rec_length = [71];

        nirs_file_number_matlab_1 = [1 ;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:71)';
            8;
            9;
            0];
     
     case 45
        N_seperate_recordings = 1;
        rec_length = [75];
        nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:75)';
            8;
            9;
            0];

        


            case 47
    N_seperate_recordings = 3;
        rec_length = [24 7 78];

        nirs_file_number_matlab_1 = [1 ;
            (10:1:19)';
            2;
            (20:1:24)';
            3;
            4;
            5;
            6;
            7;
            8;
            9;
            0];

            nirs_file_number_matlab_2 = [1 ;
            2;
            3;
            4;
            5;
            6;
            7;
            0];

            nirs_file_number_matlab_3 = [1 ;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:78)';
            8;
            9;
            0];

        
        case 48
        N_seperate_recordings = 1;
        rec_length = [77];
        nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:77)';
            8;
            9;
            0];


        case 54
    N_seperate_recordings = 1;
        rec_length = [7];

        nirs_file_number_matlab_1 = [1;    
            2;            
            3;       
            4;        
            5;         
            6;
            7;
            0];

            case 55
    N_seperate_recordings = 2;
        rec_length = [94 52];

        nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:79)';
            8;
            (80:1:89)';
            9;
            (90:1:94)';
            0];

        nirs_file_number_matlab_2 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:52)';
            6;
            7;
            8;
            9;
            0];

        case 56
    
        N_seperate_recordings = 1;
        rec_length = [140];

        nirs_file_number_matlab_1 = [1 ;
            10;
            (100:1:109)';
            11;
            (110:1:119)';
            12;
            (120:1:129)';
            13;
            (130:1:139)';
            14;
            140;
            15;
            (16:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:79)';
            8;
            (80:1:89)';
            9;
            (90:1:99)';
            0];

        case 58
    
        N_seperate_recordings = 1;
        rec_length = [111];

        nirs_file_number_matlab_1 = [1 ;
            10;
            (100:1:109)';
            11;
            (110:1:111)';
            12;
            13;
            14;
            15;
            (16:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:79)';
            8;
            (80:1:89)';
            9;
            (90:1:99)';
            0];


         case 59
            N_seperate_recordings = 1;
            rec_length = [84];

            nirs_file_number_matlab_1 = [1;
            (10:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:79)';
            8;
            (80:1:84)';
            9;
            0];



        case 60
    
        N_seperate_recordings = 1;
        rec_length = [144];

        nirs_file_number_matlab_1 = [1 ;
            10;
            (100:1:109)';
            11;
            (110:1:119)';
            12;
            (120:1:129)';
            13;
            (130:1:139)';
            14;
            (140:1:144)';
            15;
            (16:1:19)';
            2;
            (20:1:29)';
            3;
            (30:1:39)';
            4;
            (40:1:49)';
            5;
            (50:1:59)';
            6;
            (60:1:69)';
            7;
            (70:1:79)';
            8;
            (80:1:89)';
            9;
            (90:1:99)';
            0];




end %switch subjectN

%% combine nirs file number matlab together

switch N_seperate_recordings
    case 1
        nirs_file_number_matlab_1 = [nirs_file_number_matlab_1 ones(size(nirs_file_number_matlab_1,1),1) ];
       nirs_file_number_matlab = [nirs_file_number_matlab_1];

    case 2
        nirs_file_number_matlab_1 = [nirs_file_number_matlab_1 ones(size(nirs_file_number_matlab_1,1),1) ];
       nirs_file_number_matlab_2 = [nirs_file_number_matlab_2 ones(size(nirs_file_number_matlab_2,1),1)*2      ];
       nirs_file_number_matlab = [nirs_file_number_matlab_1;nirs_file_number_matlab_2];

    case 3
    nirs_file_number_matlab_1 = [nirs_file_number_matlab_1 ones(size(nirs_file_number_matlab_1,1),1) ];
       nirs_file_number_matlab_2 = [nirs_file_number_matlab_2 ones(size(nirs_file_number_matlab_2,1),1)*2      ];
       nirs_file_number_matlab_3 = [nirs_file_number_matlab_3 ones(size(nirs_file_number_matlab_3,1),1)*3      ];
       nirs_file_number_matlab = [nirs_file_number_matlab_1;nirs_file_number_matlab_2;nirs_file_number_matlab_3];

    case 4
       nirs_file_number_matlab_1 = [nirs_file_number_matlab_1 ones(size(nirs_file_number_matlab_1,1),1) ];
       nirs_file_number_matlab_2 = [nirs_file_number_matlab_2 ones(size(nirs_file_number_matlab_2,1),1)*2      ];
       nirs_file_number_matlab_3 = [nirs_file_number_matlab_3 ones(size(nirs_file_number_matlab_3,1),1)*3      ];
       nirs_file_number_matlab_4 = [nirs_file_number_matlab_4 ones(size(nirs_file_number_matlab_4,1),1)*4      ];
       nirs_file_number_matlab = [nirs_file_number_matlab_1;nirs_file_number_matlab_2;nirs_file_number_matlab_3;nirs_file_number_matlab_4];
    case 5

       nirs_file_number_matlab_1 = [nirs_file_number_matlab_1 ones(size(nirs_file_number_matlab_1,1),1) ];
       nirs_file_number_matlab_2 = [nirs_file_number_matlab_2 ones(size(nirs_file_number_matlab_2,1),1)*2      ];
       nirs_file_number_matlab_3 = [nirs_file_number_matlab_3 ones(size(nirs_file_number_matlab_3,1),1)*3      ];
       nirs_file_number_matlab_4 = [nirs_file_number_matlab_4 ones(size(nirs_file_number_matlab_4,1),1)*4      ];
       nirs_file_number_matlab_5 = [nirs_file_number_matlab_5 ones(size(nirs_file_number_matlab_5,1),1)*5      ];

       nirs_file_number_matlab = [nirs_file_number_matlab_1;nirs_file_number_matlab_2;nirs_file_number_matlab_3;nirs_file_number_matlab_4;nirs_file_number_matlab_5];


end

%% create nirs file number
nirs_file_number= ones(size(nirs_file_number_matlab,1),2);

nirs_file_number(1:rec_length(1),1) = [1:1:rec_length(1)]';
nirs_file_number(rec_length(1)+1,1) = 0;

%if there are more than 1 seperate recordings
if size(rec_length,2) >1
    for i=2:size(rec_length,2)
        %sum(rec_length(1:i-1))+2   sum(rec_length(1:i))+2 
        nirs_file_number(sum(rec_length(1:i-1))+i :sum(rec_length(1:i))+i-1 ,1) = [1:1:rec_length(i)]'; %setting 1:1:N
        nirs_file_number(sum(rec_length(1:i))+i,1) = 0; %putting zero to reflect nirs file with no number at end
        nirs_file_number(sum(rec_length(1:i-1))+i :sum(rec_length(1:i))+i,2) = i; %2nd col meas. index (1,2,3...N)
        %%%%%%%
        %nirs_file_number(rec_length(i-1)+2:rec_length(i)+rec_length(i-1)+1,1) = [1:1:rec_length(i)]';
        %nirs_file_number(rec_length(i)+rec_length(i-1)+1,1) = 0;

        %nirs_file_number(rec_length(i-1)+2:rec_length(i)+rec_length(i-1)+1,2) = i;


    end
end
%%

%test below for manual comparison
%%%%%
% nirs_file_number_1 = [1:1:18]';
% nirs_file_number_1 = [nirs_file_number_1 ; 0];
% nirs_file_number_1 = [nirs_file_number_1 ones(size(nirs_file_number_1,1),1)];
% 
% nirs_file_number_2 = [1:1:52]';
% nirs_file_number_2 = [nirs_file_number_2 ; 0];
% nirs_file_number_2 = [nirs_file_number_2 ones(size(nirs_file_number_2,1),1)*2];
% 
% nirs_file_number_3 = [1:1:18]';
% nirs_file_number_3 = [nirs_file_number_3 ; 0];
% nirs_file_number_3 = [nirs_file_number_3 ones(size(nirs_file_number_3,1),1)*3];
% 
% nirs_file_number_4 = [1:1:26]';
% nirs_file_number_4 = [nirs_file_number_4 ; 0];
% nirs_file_number_4 = [nirs_file_number_4 ones(size(nirs_file_number_4,1),1)*4];
% 
% nirs_file_number = [nirs_file_number_1;nirs_file_number_2;nirs_file_number_3;nirs_file_number_4];





end %function