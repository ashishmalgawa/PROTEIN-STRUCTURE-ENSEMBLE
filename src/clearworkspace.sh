while true; do
    read -p "This will delete output,scwrloutput,ensembles and log files.Are you Sure? (y/n)" yn
    case $yn in
        [Yy]* ) echo "Clearing workspace..."; 
        rm -r output;
        echo "output deleted...";
        rm -r scwrloutput;
        echo "scwrloutput deleted...";
        rm -r ensembles;
        echo "ensembles deleted...";
        rm *.backrub.log
        echo "log files deleted"
        break;;
        * ) 
        echo "Exiting now..";
        exit;;
    esac
done

