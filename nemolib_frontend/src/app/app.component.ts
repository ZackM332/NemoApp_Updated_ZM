import {Component, OnInit, ViewChild} from '@angular/core';
import {FormBuilder, FormControl, FormGroup, FormGroupDirective, NgForm, Validators} from '@angular/forms';
import {ErrorStateMatcher, MatPaginator, MatTableDataSource} from '@angular/material';
import {FileValidator} from 'ngx-material-file-input';
import {HttpEventType, HttpResponse} from '@angular/common/http';
import {ApiService} from './services/api.service';
import {NetworkMotifResults} from './model/NetworkMotifResults';
import {NemoProfile} from './model/NemoProfile';
import * as uuid from 'uuid';
import {FileResults} from './model/FileResults';
import * as fileSaver from 'file-saver';


/** Error when invalid control is dirty, touched, or submitted. */
export class MyErrorStateMatcher implements ErrorStateMatcher {
    isErrorState(control: FormControl | null, form: FormGroupDirective | NgForm | null): boolean {
        const isSubmitted = form && form.submitted;
        return !!(control && control.invalid && (control.dirty || control.touched || isSubmitted));
    }
}

export interface DownloadFiles {
    name: string;
    url: string;
}

@Component({
    selector: 'app-root',
    templateUrl: './app.component.html',
    styleUrls: ['./app.component.css']
})
export class AppComponent implements OnInit {
    constructor(private fb: FormBuilder, private apiService: ApiService) {
        this.mSControl.setValue('3');
        this.rSControl.setValue('10');
    }

    @ViewChild(MatPaginator) paginator: MatPaginator;

    dataSource: MatTableDataSource<NemoProfile>;
    // displayedColumns: string[] = ['label', 'nodeid', 'frequency'];
    displayedColumns: string[] = ['name', 'url'];
    downloadFileData: DownloadFiles[] = [];
    progress: { percentage: number } = {percentage: 0};
    results = '';
    maxSize = 104857600; // maximum size = 100mb
    matcher = new MyErrorStateMatcher();

    formDoc: FormGroup;

    mSControl = new FormControl('', [
        Validators.required,
        Validators.min(3),
        Validators.max(8)
    ]);
    rSControl = new FormControl('', [
        Validators.required,
        Validators.min(10)
    ]);

    probSel = '1';
    opSel = '1';
    prob = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    response: NetworkMotifResults;
    fileResponse: FileResults;
    submitted: boolean;
    resultsGet: boolean;
    currentFileUpload: boolean;
    npFinshed = false;
    graphType = '0';
    fileSystemName: string;

    ngOnInit() {
        this.formDoc = this.fb.group({
            reFile: [
                undefined, [Validators.required, FileValidator.maxContentSize(this.maxSize)]
            ]
        });
    }

    validation(): boolean {
        const motifSize = this.mSControl.value;
        if (this.probSel == null || this.opSel == null || !this.mSControl.valid || !this.rSControl.valid || !this.formDoc.valid) {
            return false;
        } else {
            for (let i = 0; i < motifSize; i++) {
                if (this.prob[i] == null || this.prob[i] <= 0 || this.prob[i] > 1) {
                    return false;
                }
            }
        }
        this.currentFileUpload = true;
        return true;
    }

    postData(): void {
        this.submitted = true;
        this.resultsGet = false;
        this.dataSource = null;
        if (this.validation()) {
            const formdata: FormData = new FormData();
            formdata.append('motifSize', this.mSControl.value);
            formdata.append('randSize', this.rSControl.value);
            formdata.append('directed', this.graphType);
            formdata.append('file', this.formDoc.get('reFile').value.files[0]);
            if (this.probSel === '1') {
                for (let i = 0; i < this.mSControl.value; i++) {
                    formdata.append('prob[]', '1.0');
                }
            } else if (this.probSel === '2') {
                for (let i = 0; i < this.mSControl.value; i++) {
                    formdata.append('prob[]', this.prob[i].toString());
                }
            }
            this.results += 'Uploading files...\n';
            this.submit(formdata);
        } else {
            this.submitted = false;
            this.currentFileUpload = false;
            alert('Please enter all parameters');
            return;
        }
    }

    submit(formdata: FormData) {
        if (this.opSel === '1') {
            this.apiService.submitNetworkMotif(formdata).subscribe(
                res => {
                    if (res.type === HttpEventType.UploadProgress) {
                        this.progress.percentage = Math.round(100 * res.loaded / res.total);
                        if (this.progress.percentage === 100) {
                            this.currentFileUpload = false;
                            this.results += 'Processing...\n';
                        }
                    } else if (res instanceof HttpResponse) {
                        this.results += 'File is completely uploaded!\n';
                        this.response = JSON.parse(res.body.toString());
                        this.results += this.response.message;
                        this.results += this.response.results;
                        this.results += '\n';
                        this.results += '-------------------------------------End of Results------------------------------------\n\n';
                        this.currentFileUpload = false;
                        this.resultsGet = true;
                        this.submitted = false;
                    }
                },
                err => {
                    alert('Sorry, an error occurred, please check your input file or contact us');
                    window.location.reload();
                }
            );
        } else if (this.opSel === '2') {
            let ids = this.formDoc.get('reFile').value.files[0].name + '_' + uuid.v4();
            ids = ids.replace(/-/g, '');
            formdata.append('uuid', ids);
            this.apiService.submitNemoProfile(formdata).subscribe(
                res => {
                    if (res.type === HttpEventType.UploadProgress) {
                        this.progress.percentage = Math.round(100 * res.loaded / res.total);
                        if (this.progress.percentage === 100) {
                            this.currentFileUpload = false;
                            this.results += 'Processing...\n';
                        }
                    } else if (res instanceof HttpResponse) {
                        this.results += 'File is completely uploaded!\n';
                        this.fileResponse = JSON.parse(res.body.toString());
                        this.results += this.fileResponse.message;
                        this.results += this.fileResponse.results;
                        this.showFileDownload();
                        this.results += '\n';
                        this.results += '-------------------------------------End of Results------------------------------------\n\n';
                        this.currentFileUpload = false;
                        this.resultsGet = true;
                        this.submitted = false;
                    }
                },
                err => {
                    alert('Sorry, an error occurred, please check your input file or contact us');
                    window.location.reload();
                }
            );
        } else if (this.opSel === '3') {
            let ids = this.formDoc.get('reFile').value.files[0].name + '_' + uuid.v4();
            ids = ids.replace(/-/g, '');
            formdata.append('uuid', ids);
            this.apiService.submitNemoCollection(formdata).subscribe(
                res => {
                    if (res.type === HttpEventType.UploadProgress) {
                        this.progress.percentage = Math.round(100 * res.loaded / res.total);
                        if (this.progress.percentage === 100) {
                            this.currentFileUpload = false;
                            this.results += 'Processing...\n';
                        }
                    } else if (res instanceof HttpResponse) {
                        this.results += 'File is completely uploaded!\n';
                        this.fileResponse = JSON.parse(res.body.toString());
                        this.results += this.fileResponse.message;
                        this.results += this.fileResponse.results;
                        this.showFileDownload();
                        this.results += '\n';
                        this.results += '-------------------------------------End of Results------------------------------------\n\n';
                        this.currentFileUpload = false;
                        this.resultsGet = true;
                        this.submitted = false;
                    }
                },
                err => {
                    alert('Sorry, an error occurred, please check your input file or contact us');
                    window.location.reload();
                }
            );
        }
    }

    downloadFile(fileName: string, url: string) {
        this.apiService.downloadFile(url).subscribe(response => {
            this.saveFile(response.body, fileName);
        });
    }

    saveFile(data: any, filename?: string) {
        const blob = new Blob([data], {type: 'text/plain'});
        fileSaver.saveAs(blob, filename);
    }

    showFileDownload() {
        if (this.fileResponse.optional == null || this.fileResponse.optional === '') {
            this.downloadFileData = this.downloadFileData.concat([{
                name: this.fileResponse.filename,
                url: this.fileResponse.url
            }]);
        } else {
            this.results += 'No File generated\n';
        }
    }

    showNemo() {
        if (this.response.optional == null || this.response.optional === '') {
            this.results += 'No NemoProfile found\n';
        } else {
            // this.results += 'Nemo profile label\n' + this.response.optional + '\n';
            const temp = this.response.optional.split('\n');
            const freqArr = [];
            for (let i = 0; i < temp.length / 2; i++) {
                const temp1 = temp[i * 2];
                const temp2 = temp[(i * 2 + 1)];
                if (temp1 != null && temp2 != null) {
                    const idArr = temp2.split('[');
                    for (let j = 1; j < idArr.length; j++) {
                        const id = idArr[j].substring(0, idArr[j].indexOf(','));
                        const fq = idArr[j].substring(idArr[j].indexOf(',') + 1, idArr[j].indexOf(']'));
                        freqArr.push({label: temp1, nodeID: id, frequency: fq});
                    }
                    idArr.shift();
                }
            }
            this.dataSource = new MatTableDataSource<NemoProfile>(freqArr);
            this.dataSource.paginator = this.paginator;
            this.npFinshed = true;
        }
    }

    applyFilter(filterValue: string) {
        this.dataSource.filter = filterValue.trim().toLowerCase();
    }

    cleanResults() {
        this.results = '';
        this.downloadFileData = [];
    }
}


