<?php

function fileUploaded($param){
    if(empty($_FILES)) {
        return FALSE;       
    }
 
    $upfile = $_FILES[$param];
    if(!file_exists($upfile['tmp_name']) || !is_uploaded_file($upfile['tmp_name'])){
        return FALSE;
    }   
    return TRUE;
}
