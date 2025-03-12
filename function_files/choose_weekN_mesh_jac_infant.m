function weekN_J = choose_weekN_mesh_jac_infant(subjectN)


switch subjectN
    case 7
        weekN_J = 31; %right
        weekN_J = 30;
    case 8
         weekN_J = 30; %right
    case 10
        weekN_J = 31; %right
        weekN_J = 30;
    case 13
         weekN_J = 32; %right
         weekN_J = 30;
    case 15
         weekN_J = 32; %right
         weekN_J = 30;
    case 20
         weekN_J = 30; %right
    case 25
         weekN_J = 31; %right
         weekN_J = 30;
    case 28
        weekN_J = 32; %right
        weekN_J = 30;
    case 33
         weekN_J = 31;%right
         weekN_J = 30;
    case 34
        weekN_J = 31; %right
        weekN_J = 30;
    case 39
         weekN_J = 31;%right
         weekN_J = 30;
    case 41
         weekN_J = 31;%right
         weekN_J = 30;
    case 44
         weekN_J = 29;%right
         weekN_J = 30;
    case 45
         weekN_J = 30;%unknown
         weekN_J = 30;
    case 48
         weekN_J = 30;%unknown
         weekN_J = 30;
    case 49
        weekN_J = 32; %right
        weekN_J = 30;

    case 55
        weekN_J = 30; %right
        weekN_J = 30;
    case 56
        weekN_J = 28; %right
        weekN_J = 30;
    case 58
        weekN_J = 30; %right
        weekN_J = 30;
    case 59
        weekN_J = 29; %right
        weekN_J = 30;
    case 60
        weekN_J = 30; %right
        weekN_J = 30;

       
end


end