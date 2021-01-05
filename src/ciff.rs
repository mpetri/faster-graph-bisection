use log::*;

pub mod proto;

pub enum CiffRecord {
    PostingsList(proto::PostingsList),
    Document {
        doc_id: u32,
        external_id: String,
        length: u32,
    },
}

pub struct Reader<'a> {
    input: protobuf::CodedInputStream<'a>,
    average_doc_len: f64,
    num_records: usize,
    num_docs: usize,
    num_postings_lists: usize,
    postings_left: usize,
    docs_left: usize,
    header: proto::Header,
}

impl<'a> Reader<'a> {
    pub fn new<T: 'a + std::io::BufRead>(input: &'a mut T) -> anyhow::Result<Reader<'a>> {
        let mut input = protobuf::CodedInputStream::from_buffered_reader(input);
        let header = input.read_message::<proto::Header>()?;
        Ok(Reader {
            input,
            average_doc_len: header.get_average_doclength(),
            postings_left: header.get_num_postings_lists() as usize,
            docs_left: header.get_num_docs() as usize,
            num_docs: header.get_num_docs() as usize,
            num_postings_lists: header.get_num_postings_lists() as usize,
            num_records: header.get_num_postings_lists() as usize + header.get_num_docs() as usize,
            header,
        })
    }

    pub fn header(&self) -> &proto::Header {
        &self.header
    }

    pub fn num_docs(&self) -> usize {
        self.num_docs
    }

    pub fn num_postings_lists(&self) -> usize {
        self.num_postings_lists
    }

    pub fn avg_doc_len(&self) -> f64 {
        self.average_doc_len
    }
}

impl<'a> ExactSizeIterator for Reader<'a> {
    fn len(&self) -> usize {
        self.num_records
    }
}

impl<'a> Iterator for Reader<'a> {
    type Item = CiffRecord;

    fn next(&mut self) -> Option<CiffRecord> {
        if self.postings_left != 0 {
            self.postings_left -= 1;
            return match self.input.read_message::<proto::PostingsList>() {
                Ok(record) => Some(CiffRecord::PostingsList(record)),
                Err(e) => {
                    error!("Error parsing CIFF postingslist: {}", e);
                    None
                }
            };
        }
        if self.docs_left != 0 {
            self.docs_left -= 1;
            return match self.input.read_message::<proto::DocRecord>() {
                Ok(record) => {
                    Some(CiffRecord::Document {
                        doc_id: record.get_docid() as u32, // todo fix this...
                        external_id: record.get_collection_docid().to_string(),
                        length: record.get_doclength() as u32, // todo fix this...
                    })
                }
                Err(e) => {
                    error!("Error parsing CIFF document: {}", e);
                    None
                }
            };
        }
        None
    }
}


pub struct Writer<'a> {
    output: protobuf::CodedOutputStream<'a>,
}

impl<'a> Writer<'a> {
    pub fn from_reader<T: 'a + std::io::Write>(output: &'a mut T,reader: &Reader<'a>,description: &str) -> Result<Writer<'a>, protobuf::ProtobufError> {
        let header = reader.header();
        let description = format!("{} {}",header.get_description(),description);
        Self::new(
            output,
            header.get_num_postings_lists(),
            header.get_total_postings_lists(),
            header.get_total_terms_in_collection(),
            header.get_num_docs(),
            header.get_total_docs(),
            header.get_average_doclength(),
            &description
        )
    }

    pub fn new<T: 'a + std::io::Write>(
        output: &'a mut T,
        num_postings_lists: i32,
        total_postings_lists: i32,
        total_terms_in_collection: i64,
        num_docs: i32,
        total_docs: i32,
        avg_doc_len: f64,
        description: &str,
    ) -> Result<Writer<'a>, protobuf::ProtobufError> {
        let mut writer = Writer {
            output: protobuf::CodedOutputStream::new(output),
        };

        // write the header
        let mut header = proto::Header::default();
        header.set_version(1);
        header.set_description(description.into());
        header.set_num_postings_lists(num_postings_lists);
        header.set_total_postings_lists(total_postings_lists);
        header.set_total_terms_in_collection(total_terms_in_collection);
        header.set_num_docs(num_docs);
        header.set_total_docs(total_docs);
        header.set_average_doclength(avg_doc_len);

        writer.output.write_message_no_tag(&header)?;

        Ok(writer)
    }

    pub fn write_document(
        &mut self,
        docid: i32,
        external_id: String,
        len: i32,
    ) -> Result<(), protobuf::ProtobufError> {
        let mut document = proto::DocRecord::default();
        document.set_docid(docid);
        document.set_collection_docid(external_id);
        document.set_doclength(len);
        self.output.write_message_no_tag(&document)
    }

    pub fn write_postingslist(
        &mut self,
        postings_list: proto::PostingsList,
    ) -> Result<(), protobuf::ProtobufError> {
        self.output.write_message_no_tag(&postings_list)
    }

    pub fn flush(&mut self) -> Result<(), protobuf::ProtobufError> {
        self.output.flush()
    }
}